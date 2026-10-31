# Performance checks and the regrid parity check.
#   julia --project=prototype/env prototype/bench.jl
using GridInterface, BenchmarkTools, Test, SparseArrays, LinearAlgebra
import GridInterface as GrI
using DiscreteGlobalGrids
import DiscreteGlobalGrids as DGG
using GridInterface: cell_extent, cellat   # DGG exports its own (a lon/lat box, a typed id)
import ConservativeRegridding as CR                       # regrids GridInterface grids directly (Phase 7)
using GridInterface.Structured
import Extents, GeometryOps as GO, GeometryOpsCore as GOC
const FB = DGG.Fallbacks

fmt(t) = t < 1e-6 ? string(round(t * 1e9; digits=1), " ns") : t < 1e-3 ? string(round(t * 1e6; digits=2), " µs") :
         t < 1 ? string(round(t * 1e3; digits=2), " ms") : string(round(t; digits=2), " s")
section(s) = (println("\n=== ", s, " ==="); flush(stdout))

sweep_centroid(g) = (s = 0.0; for i in 1:ncells(g); s += cell_centroid(g, i)[1]; end; s)
sweep_extent(g) = (s = 0.0; for i in 1:ncells(g); s += cell_extent(g, i).radius; end; s)
sweep_identity(g) = (s = 0; for i in 1:ncells(g); s += localindex(g, cellindex(g, i))::Int; end; s)
direct_centroid(g) = (s = 0.0; for i in 1:DGG.ncells(g); s += DGG.cell_centroid(g, DGG.cellindex(g, i))[1]; end; s)
direct_set_centroid(set) = (s = 0.0; sys = set.system; for c in set.cells; s += DGG.cell_centroid(DGG.levelgrid(sys, DGG.level(c)), c)[1]; end; s)

"CR's standard global lon/lat test grid: the corner matrix (or the axes) as a `Structured` grid with its seams declared."
function lonlat_grid(nlon, nlat; manifold=GOC.Spherical(), cellbased=true)
    lons = range(-180, 180; length=nlon + 1); lats = range(-90, 90; length=nlat + 1)
    seams = GrI.iscurved(manifold) ? (; periodic=(true, false), collapsed=(true, true)) : (; periodic=(false, false), collapsed=(false, false))
    cellbased || return Structured.RegularGrid(manifold, collect(lons), collect(lats); seams...)
    return Structured.CellBasedGrid(manifold, [GO.UnitSphereFromGeographic()((lo, la)) for lo in lons, la in lats]; seams...)
end
"The same corners as a bare matrix: what `CR.Regridder` turns into a `CellBasedGrid` itself."
lonlat_points(nlon, nlat) = [GO.UnitSphereFromGeographic()((lo, la)) for lo in range(-180, 180; length=nlon + 1), la in range(-90, 90; length=nlat + 1)]
# a grid paired with GridInterface's generic R-tree over its cell extents, instead of its own cursor
rtree(g) = GridAndTree(g, GrI.fallback_treeify(GO.best_manifold(g), g, 16))

section("1. per-cell geometry: level grids and a multi-level set (type stability, allocation)")
hp = DGG.HEALPixSystem()
set = DGG.query(hp, DGG.MultiOrderCoverage(Extents.Extent(X=(-10.0, 30.0), Y=(35.0, 60.0))); level=10)
mls = set                                  # a MultiOrderCellSet is a GridInterface multi-level grid
println("HEALPix multi-order set: ", ncells(mls), " members, levels ", extrema(GrI.level.(mls.cells)), ", antichain ", GrI.check_antichain(mls))
for (name, g, direct) in (("HEALPix level 6", GrI.levelgrid(hp, 6), direct_centroid),
                          ("IGeo7 level 5", GrI.levelgrid(DGG.IGeo7System(), 5), direct_centroid),
                          ("H3 level 3", GrI.levelgrid(DGG.H3System(), 3), direct_centroid),
                          ("HEALPix multi-level set", mls, nothing))
    @inferred cell_centroid(g, 1); @inferred cell_extent(g, 1); @inferred cellindex(g, 1)
    n = ncells(g)
    tc = @belapsed sweep_centroid($g); ac = @ballocated sweep_centroid($g)
    te = @belapsed sweep_extent($g); ae = @ballocated sweep_extent($g)
    ti = @belapsed sweep_identity($g); ai = @ballocated sweep_identity($g)
    td = direct === nothing ? (@belapsed direct_set_centroid($set)) : (@belapsed $direct($g))
    println(rpad(name, 26), "n=", n, "  @inferred ok  centroid sweep ", fmt(tc / n), "/cell ", ac, " B (direct DGG ", fmt(td / n), "/cell)",
        "  extent ", fmt(te / n), "/cell ", ae, " B  identity ", fmt(ti / n), "/cell ", ai, " B")
end
println("levelgrid return types (HEALPix, IGeo7, H3, A5): ",
    [Base.return_types(GrI.levelgrid, Tuple{typeof(s),Int}) for s in (hp, DGG.IGeo7System(), DGG.H3System(), DGG.A5System())])

section("2. treeify on a level grid is O(1)")
for l in (3, 6, 9, 12)
    g = GrI.levelgrid(hp, l)
    t = @belapsed treeify($g); a = @ballocated treeify($g)
    println("HEALPix level $l (", ncells(g), " cells): treeify ", fmt(t), ", ", a, " B")
end
tm = @belapsed treeify($mls); am = @ballocated treeify($mls)
println("multi-level set: treeify ", fmt(tm), ", ", am, " B")
let g = GrI.levelgrid(DGG.A5System(), 4)
    t = @belapsed treeify($g); a = @ballocated treeify($g)
    println("A5 level 4 (selection mode, ", ncells(g), " cells): treeify ", fmt(t), ", ", a, " B  (O(n): materialises the selection)")
end

let cg = GrI.levelgrid(DGG.CopernicusDEMSystem(90), 1)          # 6.9e10 cells
    mt = treeify(cg, GrI.ForQueries())
    tiny = GO.UnitSpherical.SphericalCap(GO.UnitSphereFromGeographic()((10.123, 45.456)), 1e-7)
    hits(t, q) = (n = Ref(0); GO.SpatialTreeInterface.depth_first_search(e -> Extents.intersects(e, q), t) do i; n[] += 1; nothing; end; n[])
    nl = hits(mt, tiny)
    t = @belapsed $hits($mt, $tiny)
    println("CopernicusDEM(90) L1: treeify(ForQueries()) is ", nameof(typeof(mt)), "; a point query reaches ", nl, " leaves in ", fmt(t))
end

section("3. cellat on analytic grids")
for (name, g) in (("HEALPix 8", GrI.levelgrid(hp, 8)), ("IGeo7 6", GrI.levelgrid(DGG.IGeo7System(), 6)),
                  ("H3 5", GrI.levelgrid(DGG.H3System(), 5)), ("S2 8", GrI.levelgrid(DGG.S2System(), 8)),
                  ("ISEA4R 6", GrI.levelgrid(DGG.ISEA4RSystem(), 6)), ("A5 5", GrI.levelgrid(DGG.A5System(), 5)),
                  ("HEALPix multi-level set", mls))
    p = cell_centroid(g, ncells(g) ÷ 3)
    @assert cellat(g, p) == ncells(g) ÷ 3
    t = @belapsed cellat($g, $p); a = @ballocated cellat($g, $p)
    println(rpad(name, 26), "analytic=", location_is_analytic(g), "  cellat ", fmt(t), ", ", a, " B")
end
let g = mls, p = cell_centroid(mls, 100), gt = GridAndTree(mls, GrI.fallback_treeify(GO.Spherical(radius=1.0), mls, 16))
    GrI.fallback_cellat(GO.Spherical(radius=1.0), g, p)
    t = @belapsed GrI.fallback_cellat(GO.Spherical(radius=1.0), $g, $p)
    t2 = @belapsed cellat($gt, $p)
    println("multi-level set fallback cellat via hierarchy cursor ", fmt(t), "; via cached R-tree (GridAndTree) ", fmt(t2))
end

section("4. regrid parity: DGG's own cursors through CR's Regridder vs the generic R-tree")
function parity(name, dst, src)
    R_dgg = CR.Regridder(dst, src; threaded=false)                      # DGG's HierarchicalGridCursor (the provider's treeify)
    R_gi = CR.Regridder(rtree(dst), rtree(src); threaded=false)         # GridInterface's generic R-tree
    t_dgg = @elapsed CR.Regridder(dst, src; threaded=false)
    t_gi = @elapsed CR.Regridder(rtree(dst), rtree(src); threaded=false)
    A, B = R_dgg.intersections, R_gi.intersections
    println(rpad(name, 34), size(A), " nnz ", nnz(A), "  identical: ", A == B, "  max|Δ| ", maximum(abs.(A - B); init=0.0),
        "  areas identical: ", R_dgg.dst_areas == R_gi.dst_areas && R_dgg.src_areas == R_gi.src_areas,
        "  time cursor ", fmt(t_dgg), " / R-tree ", fmt(t_gi))
    return A == B
end
parity("HEALPix 5 ← IGeo7 4", GrI.levelgrid(hp, 5), GrI.levelgrid(DGG.IGeo7System(), 4))
parity("S2 5 ← H3 2", GrI.levelgrid(DGG.S2System(), 5), GrI.levelgrid(DGG.H3System(), 2))
parity("ISEA4R 4 ← A5 3 (selection mode)", GrI.levelgrid(DGG.ISEA4RSystem(), 4), GrI.levelgrid(DGG.A5System(), 3))
pg = DGG.subtree(hp, DGG.cellindex(GrI.levelgrid(hp, 1), 7), 6)
parity("IGeo7 4 ← HEALPix partial grid", GrI.levelgrid(DGG.IGeo7System(), 4), pg)

section("5. multi-level set straight into CR's Regridder")
dst = GrI.levelgrid(DGG.IGeo7System(), 5)
R_set = CR.Regridder(dst, mls; threaded=false)
t_set = @elapsed CR.Regridder(dst, mls; threaded=false)
rt = GrI.fallback_treeify(GO.Spherical(radius=1.0), mls, 16)
R_rt = CR.Regridder(dst, GridAndTree(mls, rt); threaded=false)
t_rt = @elapsed CR.Regridder(dst, GridAndTree(mls, rt); threaded=false)
println("IGeo7 5 ← HEALPix multi-level set: ", size(R_set.intersections), " cursor ", fmt(t_set), "; R-tree ", fmt(t_rt),
    "; identical to R-tree: ", R_set.intersections == R_rt.intersections)
expanded = DGG.PartialGrid(DGG.CellVector(set))
t_exp = @elapsed (R_exp = CR.Regridder(dst, expanded; threaded=false))
rows, cols = GrI.aggregation_pairs(mls)
lg_ref = GrI.levelgrid(hp, mls.reference_level)
# expanded columns are the partial grid's own positions; map reference positions to them
colpos = [localindex(expanded, DGG.cellindex(lg_ref, r)) for r in rows]
P = sparse(colpos, cols, 1.0, ncells(expanded), ncells(mls))
gap = maximum(abs.(R_exp.intersections * P - R_set.intersections)) / maximum(R_set.intersections)
println("today's expansion path: ", size(R_exp.intersections), " in ", fmt(t_exp), "; R_expanded·P vs R_set max rel gap ", gap)

section("6. regrid parity: CR's bare inputs vs the Structured grids they become (bit-identical)")
function cr_parity(name, dst_bare, src_bare, dst_grid, src_grid)
    R_cr = CR.Regridder(dst_bare, src_bare; threaded=false)                 # a corner matrix / axis pair: `asgrid` builds the grid
    R_gi = CR.Regridder(dst_grid, src_grid; threaded=false)                 # the Structured grids built explicitly
    t_cr = @elapsed CR.Regridder(dst_bare, src_bare; threaded=false)
    t_gi = @elapsed CR.Regridder(dst_grid, src_grid; threaded=false)
    A, B = R_cr.intersections, R_gi.intersections
    same = A == B && R_cr.dst_areas == R_gi.dst_areas && R_cr.src_areas == R_gi.src_areas
    println(rpad(name, 50), size(A), " nnz ", nnz(A), "  bit-identical: ", same, "  max|Δ| ", maximum(abs.(A - B); init=0.0),
        "  time bare ", fmt(t_cr), " / grids ", fmt(t_gi))
    return same
end
cr_parity("lon/lat 36×18 ← lon/lat 72×36 (CellBasedGrid)", lonlat_points(36, 18), lonlat_points(72, 36), lonlat_grid(36, 18), lonlat_grid(72, 36))
# A spherical RegularGrid has the same geodesic quads as the CellBasedGrid of the same corners.
let rg = lonlat_grid(36, 18; cellbased=false), ll72 = lonlat_grid(72, 36)
    R_cr = CR.Regridder(lonlat_grid(36, 18), ll72; threaded=false)
    R_gi = CR.Regridder(rg, ll72; threaded=false)
    println(rpad("RegularGrid lon/lat 36×18 ← lon/lat 72×36", 50), size(R_gi.intersections), " nnz ", nnz(R_gi.intersections),
        "  bit-identical to CellBasedGrid 36×18: ", R_cr.intersections == R_gi.intersections && R_cr.dst_areas == R_gi.dst_areas)
end
let xd = collect(0.0:1.0:20.0), yd = collect(0.0:1.0:10.0), xs = collect(range(0.0, 20.0; length=31)), ys = collect(range(0.0, 10.0; length=16))
    cr_parity("RegularGrid planar 20×10 ← RegularGrid planar 30×15", (xd, yd), (xs, ys),
        Structured.RegularGrid(GOC.Planar(), xd, yd), Structured.RegularGrid(GOC.Planar(), xs, ys))
end
let dst = GrI.levelgrid(hp, 3), m = GOC.Spherical(; radius=1.0)   # DGG's unit sphere
    R_cr = CR.Regridder(m, dst, lonlat_points(36, 18); threaded=false)
    R_gi = CR.Regridder(m, dst, lonlat_grid(36, 18; manifold=m); threaded=false)
    println(rpad("HEALPix L3 ← lon/lat 36×18 (CellBasedGrid)", 50), size(R_cr.intersections), " nnz ", nnz(R_cr.intersections),
        "  bit-identical: ", R_cr.intersections == R_gi.intersections && R_cr.src_areas == R_gi.src_areas)
end
