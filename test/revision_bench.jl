# Benchmarks for the round-one revisions (prototype/bench.jl is unchanged).
#   julia --project=prototype/env test/revision_bench.jl
using GridInterface, BenchmarkTools, SparseArrays
import GridInterface as GrI
const PROTO = joinpath(@__DIR__, "..", "prototype")
include(joinpath(PROTO, "adapters", "DGGAdapter.jl")); using .DGGAdapter
include(joinpath(PROTO, "adapters", "CRAdapter.jl")); using .CRAdapter
import Extents, GeometryOps as GO, GeometryOpsCore as GOC
using GeometryOps.UnitSpherical: UnitSphericalPoint
const STI = GrI.STI

fmt(t) = t < 1e-6 ? string(round(t * 1e9; digits=1), " ns") : t < 1e-3 ? string(round(t * 1e6; digits=2), " µs") :
         t < 1 ? string(round(t * 1e3; digits=2), " ms") : string(round(t; digits=2), " s")
section(s) = (println("\n=== ", s, " ==="); flush(stdout))
hp = DGG.HEALPixSystem()

section("A. analytic HEALPix vertices vs the geometric vertex table")
let g = GrI.levelgrid(hp, 3)          # warm up both paths
    GrI.fallback_vertex_table(GOC.best_manifold(g), g); GrI.vertex_numbering(g); cellat(dual(g), cell_centroid(g, 5))
    for i in 1:ncells(g), k in 1:4; vertex_cells(g, vertex_key(g, i, k)); end
end
for l in (5, 6, 8)
    g = GrI.levelgrid(hp, l)
    GrI.clear_caches!()
    t_tab = @elapsed GrI.fallback_vertex_table(GOC.best_manifold(g), g)
    tab = GrI.fallback_vertex_table(GOC.best_manifold(g), g)
    # analytic incidence per corner, no table
    t_key = @elapsed for i in 1:ncells(g), k in 1:4; vertex_cells(g, vertex_key(g, i, k)); end
    i = ncells(g) ÷ 3
    key = vertex_key(g, i, 2)
    t_vc_an = @belapsed vertex_cells($g, $key)
    v = tab.cell_verts[tab.cell_ptr[i]+1]
    t_vc_fb = @belapsed view($tab.vert_cells, $tab.vert_ptr[$v]:$tab.vert_ptr[$v+1]-1)
    t_num = @elapsed nvertices(g); nv = nvertices(g)
    t_numv = @elapsed GrI.vertex_numbering(g)
    D = dual(g); p = cell_centroid(g, i)
    cellat(D, p); t_dual = @belapsed cellat($D, $p)
    same = Set(sort(collect(vertex_cells(g, vertex_key(g, c, k)))) for c in 1:ncells(g) for k in 1:4) ==
           Set(sort(tab.vert_cells[tab.vert_ptr[w]:tab.vert_ptr[w+1]-1]) for w in 1:length(tab.points))
    println("HEALPix L$l ($(ncells(g)) cells, $nv vertices): table build ", fmt(t_tab), "; analytic keys for every corner ",
        fmt(t_key), " (", fmt(t_key / (4ncells(g))), "/corner, no table); vertex_cells analytic ", fmt(t_vc_an),
        " vs table ", fmt(t_vc_fb), "; dense numbering derived in ", fmt(t_numv), "; cellat(dual) ", fmt(t_dual),
        "; same incidence as the table: ", same)
end

section("B. Index(Global()) neighbours on a subset (round one: a linear scan)")
lg = GrI.levelgrid(hp, 8)
s = subset(lg, 1000:5000)
neighbors(s, 7, Index(Global())) |> collect
t1 = @belapsed collect(neighbors($s, 7, Index(Global())))
t2 = @belapsed collect(neighbors($s, 7, Index(Local())))
t3 = @belapsed collect(neighbors($s, 7))
println("subset of HEALPix L8: Index(Global()) ", fmt(t1), "; Index(Local()) ", fmt(t2), " (with OUTSIDE sentinel); clipped ", fmt(t3))

section("C. CopernicusDEM(90): the generic cursor over 64,800 roots")
cd = DGG.CopernicusDEMSystem(90)
cg = GrI.levelgrid(cd, 1)
t_ri = @elapsed GrI.root_index(cd)
t_tr = @belapsed GrI.hierarchy_tree($cg)
ct = GrI.hierarchy_tree(cg)
count_leaves(t, q) = (n = Ref(0); STI.depth_first_search(e -> Extents.intersects(e, q), t) do i; n[] += 1; nothing; end; n[])
p0 = GO.UnitSphereFromGeographic()((10.123, 45.456))
tiny = GO.UnitSpherical.SphericalCap(p0, 1e-7)
nl = count_leaves(ct, tiny)
t_q = @belapsed count_leaves($ct, $tiny)
mt = treeify(cg); count_leaves(mt, tiny)
t_qm = @belapsed count_leaves($mt, $tiny)
roots = GrI.rootcells(cd)
t_scan = @elapsed count(r -> Extents.intersects(GrI.node_extent(cd, r), tiny), roots)
println("ncells L1 = ", ncells(cg), "; root index built once in ", fmt(t_ri), "; hierarchy_tree ", fmt(t_tr),
    "; a point query reaches ", nl, " leaves in ", fmt(t_q), " through the root index (DGG's MemoBlockCursor: ", fmt(t_qm),
    "); without the index the root fan-out alone scans 64,800 root extents in ", fmt(t_scan))
println("treeify(cg) (provider tree wins): ", nameof(typeof(treeify(cg))), "; chunks(cg): ", length(chunks(cg)), " (one per 1° tile)")

section("D. split_weight through every wrapper and the CR bridge")
g6 = GrI.levelgrid(hp, 6)
for (name, tr) in (("DGG cursor", treeify(g6)), ("generic cursor", GrI.hierarchy_tree(g6)))
    gt = GridAndTree(g6, tr)
    println(rpad(name, 16), "split_weight: tree ", GrI.split_weight(tr), ", GridAndTree ", GrI.split_weight(gt),
        ", CRTree ", Trees.split_weight(cr_tree(gt)), " (ncells ", ncells(g6), ")")
end
import Healpix; adapt_healpix!(Healpix)
hm = Healpix.HealpixMap{Float64,Healpix.RingOrder}(16)
println("HealpixMap nside 16 via the bridge: ", Trees.split_weight(cr_tree(hm)), " (ncells ", ncells(hm), "; round one: 1)")

section("E. chunks: lazy storage partitions and cheap extents")
for l in (8, 12)
    g = GrI.levelgrid(hp, l)
    chunks(g); a = @allocated chunks(g); t = @belapsed chunks($g)
    ch = chunks(g)
    te = @belapsed GrI.grid_extent($ch[17])
    println("HEALPix L$l ($(ncells(g)) cells): chunks ", fmt(t), ", ", a, " B; ", length(ch), " chunks of ", ncells(ch[1]),
        " cells; chunk extent ", fmt(te), "; chunkat ", chunkat(ch, ncells(g) ÷ 2))
end
ll = lonlat_grid(1440, 720)
tc = chunks(ll)
te = @belapsed GrI.grid_extent($tc[100])
tt = @belapsed GrI.subtree($ll, $tc[100])
tl = @belapsed treeify($tc[100])
println("CR lon/lat 1440×720: ", length(tc), " TileIndices chunks; extent ", fmt(te), "; subtree ", fmt(tt), "; treeify(tile) ", fmt(tl),
    "; chunkat O(1): ", chunkat(tc, 500_000))

section("F. halo / border on a rooted chunk: DGG's subtree engine vs the generic merge")
g10 = GrI.levelgrid(hp, 10)
c = chunks(g10)[100]
r = parentpositions(c)
halo(c); t_dgg = @belapsed halo($c)
t_gen = @belapsed invoke(GrI.halo, Tuple{Any,AbstractVector{<:Integer},Integer}, $g10, $(collect(r)), 1)
border(c); t_bd = @belapsed border($c)
println("HEALPix L10 chunk (", ncells(c), " cells): halo ", fmt(t_dgg), " via DGG vs ", fmt(t_gen), " generic; border ", fmt(t_bd),
    "; same: ", halo(c) == invoke(GrI.halo, Tuple{Any,AbstractVector{<:Integer},Integer}, g10, collect(r), 1))

section("G. query predicates (lon/lat box read as DGG does: densified boundary, great-circle joins)")
g8 = GrI.levelgrid(hp, 8)
box = Extents.Extent(X=(-10.0, 30.0), Y=(35.0, 60.0))
for P in (Intersects, CentroidIn, Within, Contains)
    q = query(g8, P(box))
    t = @elapsed query(g8, P(box))
    println(rpad(string(P), 11), length(q), " cells in ", fmt(t))
end

section("H. the four-method grid regrids")
struct FourLL <: GrI.AbstractGrid
    nlon::Int
    nlat::Int
end
GrI.ncells(g::FourLL) = g.nlon * g.nlat
GOC.best_manifold(::FourLL) = GOC.Spherical(; radius=1.0)
GrI.cellshape(::FourLL) = GrI.ConvexCells()
function GrI.cell_boundary(g::FourLL, i::Integer)
    a = (i - 1) % g.nlon; b = (i - 1) ÷ g.nlon
    x0 = -180 + 360a / g.nlon; x1 = -180 + 360(a + 1) / g.nlon
    y0 = -90 + 180b / g.nlat; y1 = -90 + 180(b + 1) / g.nlat
    f = GO.UnitSphereFromGeographic()
    return (f((x0, y0)), f((x1, y0)), f((x1, y1)), f((x0, y1)))
end
src4 = FourLL(72, 36); src_cr = lonlat_grid(72, 36; manifold=GOC.Spherical(; radius=1.0))
dst = GrI.levelgrid(hp, 3)
R4 = regridder(dst, src4; threaded=false); Rc = CR.Regridder(dst, src_cr; threaded=false)
println("HEALPix L3 ← four-method lon/lat 72×36: nnz ", nnz(R4.intersections), ", identical to CR's own lon/lat grid: ",
    R4.intersections ≈ Rc.intersections, " (max|Δ| ", maximum(abs.(R4.intersections - Rc.intersections)), ")")
