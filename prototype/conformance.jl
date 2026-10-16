# Conformance on every adapted grid.
#   julia --project=prototype/env prototype/conformance.jl
using GridInterface, GridInterfaceConformance
import GridInterface as GrI, Extents
import GeometryOps as GO, GeometryOpsCore as GOC, GeoInterface as GI
include(joinpath(@__DIR__, "adapters", "DGGAdapter.jl")); using .DGGAdapter
include(joinpath(@__DIR__, "adapters", "CRAdapter.jl")); using .CRAdapter
import Healpix; adapt_healpix!(Healpix)

report(rep, t) = (print_report(rep); println("  (", round(t; digits=1), " s)"); flush(stdout))

# DGG systems: hierarchy laws + one level grid + a sample multi-level set
for (name, sys, l) in dgg_systems()
    t = @elapsed rep = check_system(sys; name, grid_level=l)
    report(rep, t)
end
# DGG partial grids and real multi-order coverages
for (name, sys) in (("HEALPix", DGG.HEALPixSystem()), ("IGeo7", DGG.IGeo7System()), ("H3", DGG.H3System()))
    pg = DGG.subtree(sys, GrI.rootcells(sys)[2], 4)
    t = @elapsed rep = check_grid(pg; name="$name PartialGrid (subtree of root 2 to level 4)")
    report(rep, t)
    set = DGG.query(sys, DGG.MultiOrderCoverage(Extents.Extent(X=(-10.0, 30.0), Y=(35.0, 60.0))); level=6)
    t = @elapsed rep = check_multilevel(dgg_multilevel(set); name="$name MultiOrderCoverage Europe, level ≤ 6")
    report(rep, t)
end
# CR grids, through their Structured twins (what any consumer given a CR grid builds)
warp = [GO.UnitSphereFromGeographic()((lo + 6sind(3la), la + 4sind(2lo))) for lo in range(-60, 60; length=25), la in range(-40, 50; length=19)]
polys = [GI.Polygon([[(i + 0.2j, j), (i + 1 + 0.2j, j), (i + 1.2 + 0.2j, j + 1), (i + 0.2 + 0.2j, j + 1), (i + 0.2j, j)]]) for i in 0:7, j in 0:5]
for (name, g) in (
        ("CR CellBasedGrid lon/lat 36x18 (Spherical)", lonlat_grid(36, 18)),
        ("CR CellBasedGrid warped curvilinear 24x18 (Spherical)", structured(Trees.CellBasedGrid(GOC.Spherical(), warp))),
        ("CR RegularGrid 20x10 (Planar)", structured(Trees.RegularGrid(GOC.Planar(), collect(0.0:1.0:20.0), collect(0.0:1.0:10.0)))),
        ("CR RegularGrid lon/lat 36x18 (Spherical)", lonlat_grid(36, 18; cellbased=false)),
        ("CR ExplicitPolygonGrid sheared 8x6 (Planar)", structured(Trees.ExplicitPolygonGrid(GOC.Planar(), polys))),
        ("Healpix.HealpixMap nside 4 via CR ext", Healpix.HealpixMap{Float64,Healpix.RingOrder}(4)))
    t = @elapsed rep = check_grid(g; name)
    report(rep, t)
end
# The laziness law at scale: no verb's bytes scale between a small level and a huge one
for (name, sys, ls) in (("HEALPix", DGG.HEALPixSystem(), (4, 12)), ("CopernicusDEM(90)", DGG.CopernicusDEMSystem(90), (0, 1)))
    small, large = GrI.levelgrid(sys, ls[1]), GrI.levelgrid(sys, ls[2])
    t = @elapsed rep = check_laziness(small, large; name="$name laziness L$(ls[1]) / L$(ls[2])")
    report(rep, t)
end
