# Conformance on every DGG grid (DiscreteGlobalGrids is a GridInterface
# provider itself) and every adapted CR grid.
#   julia --project=prototype/env prototype/conformance.jl
using GridInterface, GridInterfaceConformance
import GridInterface as GrI, Extents
import GeometryOps as GO, GeometryOpsCore as GOC, GeoInterface as GI
using DiscreteGlobalGrids
import DiscreteGlobalGrids as DGG
import ConservativeRegridding as CR
using GridInterface.Structured
import Healpix
const HPExt = Base.get_extension(CR, :ConservativeRegriddingHealpixExt)      # HealpixMapGrid: CR's own provider

"CR's standard global lon/lat test grid as a `Structured` grid with its seams declared (what `CR.Regridder` builds from the corners)."
function lonlat_grid(nlon, nlat; manifold=GOC.Spherical(), cellbased=true)
    lons = range(-180, 180; length=nlon + 1); lats = range(-90, 90; length=nlat + 1)
    seams = GrI.iscurved(manifold) ? (; periodic=(true, false), collapsed=(true, true)) : (; periodic=(false, false), collapsed=(false, false))
    cellbased || return Structured.RegularGrid(manifold, collect(lons), collect(lats); seams...)
    return Structured.CellBasedGrid(manifold, [GO.UnitSphereFromGeographic()((lo, la)) for lo in lons, la in lats]; seams...)
end

report(rep, t) = (print_report(rep); println("  (", round(t; digits=1), " s)"); flush(stdout))

"The registered DGG systems plus CopernicusDEM(90), with a small level for conformance runs (H3 at an odd level too)."
const DGG_SYSTEMS = [
    ("HEALPix", HEALPixSystem(), 3),
    ("IGeo7", IGeo7System(), 3),
    ("H3", H3System(), 2),
    ("H3 (odd level)", H3System(), 1),
    ("ISEA4R", ISEA4RSystem(), 3),
    ("S2", S2System(), 3),
    ("A5", A5System(), 2),
    ("CopernicusDEM(90)", CopernicusDEMSystem(90), 0),
]

# A DGG level grid is a GridInterface grid as it stands (no adapter), and DGG's
# cursors answer every `treeify` intent: the hierarchy cursor on the true
# hierarchies, CopernicusDEM's block cursor (never the generic hierarchy cursor).
let g = GrI.levelgrid(DGG.HEALPixSystem(), 4), cg = GrI.levelgrid(DGG.CopernicusDEMSystem(90), 1)
    for how in (GrI.ForQueries(), GrI.ForSweeps(), GrI.Balanced(), 16)
        @assert treeify(g, how) isa DGG.Engine.HierarchicalGridCursor
        @assert treeify(cg, how) isa DGG.CopernicusDEM.MemoBlockCursor
    end
    t = @elapsed rep = check_grid(g; name="HEALPix L4 (DGG level grid, no adapter)")
    report(rep, t)
end
# DGG systems: hierarchy laws + one level grid + a sample multi-level set
for (name, sys, l) in DGG_SYSTEMS
    t = @elapsed rep = check_system(sys; name, grid_level=l)
    report(rep, t)
end
# DGG partial grids and real multi-order coverages
for (name, sys) in (("HEALPix", DGG.HEALPixSystem()), ("IGeo7", DGG.IGeo7System()), ("H3", DGG.H3System()))
    pg = DGG.subtree(sys, GrI.rootcells(sys)[2], 4)
    t = @elapsed rep = check_grid(pg; name="$name PartialGrid (subtree of root 2 to level 4)")
    report(rep, t)
    set = DGG.query(sys, DGG.MultiOrderCoverage(Extents.Extent(X=(-10.0, 30.0), Y=(35.0, 60.0))); level=6)
    t = @elapsed rep = check_multilevel(set; name="$name MultiOrderCoverage Europe, level ≤ 6")
    report(rep, t)
end
# CR's grids (GridInterface.Structured grids since Phase 7: what `CR.Regridder` builds from bare corners, axes or polygons)
warp = [GO.UnitSphereFromGeographic()((lo + 6sind(3la), la + 4sind(2lo))) for lo in range(-60, 60; length=25), la in range(-40, 50; length=19)]
polys = [GI.Polygon([[(i + 0.2j, j), (i + 1 + 0.2j, j), (i + 1.2 + 0.2j, j + 1), (i + 0.2 + 0.2j, j + 1), (i + 0.2j, j)]]) for i in 0:7, j in 0:5]
for (name, g) in (
        ("CR CellBasedGrid lon/lat 36x18 (Spherical)", lonlat_grid(36, 18)),
        ("CR CellBasedGrid warped curvilinear 24x18 (Spherical)", CR.CellBasedGrid(GOC.Spherical(), warp)),
        ("CR RegularGrid 20x10 (Planar)", CR.RegularGrid(GOC.Planar(), collect(0.0:1.0:20.0), collect(0.0:1.0:10.0))),
        ("CR RegularGrid lon/lat 36x18 (Spherical)", lonlat_grid(36, 18; cellbased=false)),
        ("CR ExplicitPolygonGrid sheared 8x6 (Planar)", CR.ExplicitPolygonGrid(GOC.Planar(), polys)),
        ("Healpix.HealpixMap nside 4 (CR's HealpixMapGrid)", HPExt.HealpixMapGrid(Healpix.HealpixMap{Float64,Healpix.RingOrder}(4))))
    t = @elapsed rep = check_grid(g; name)
    report(rep, t)
end
# The laziness law at scale: no verb's bytes scale between a small level and a huge one
for (name, sys, ls) in (("HEALPix", DGG.HEALPixSystem(), (4, 12)), ("CopernicusDEM(90)", DGG.CopernicusDEMSystem(90), (0, 1)))
    small, large = GrI.levelgrid(sys, ls[1]), GrI.levelgrid(sys, ls[2])
    t = @elapsed rep = check_laziness(small, large; name="$name laziness L$(ls[1]) / L$(ls[2])")
    report(rep, t)
end
# The full grid laws on a 6.9e10-cell grid: every law must sample, never enumerate (about 12 s cold)
let g = GrI.levelgrid(DGG.CopernicusDEMSystem(90), 1)
    t = @elapsed rep = check_grid(g; name="CopernicusDEM(90) L1 (6.9e10 cells, sampled laws)")
    report(rep, t)
end
