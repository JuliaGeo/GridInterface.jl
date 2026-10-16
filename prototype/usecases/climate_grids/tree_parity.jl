# Tree parity: for every grid of CR's test/trees/grids.jl, the Structured
# quadtree built from the same corners yields the same node extents and leaf
# indices as CR's own TopDownQuadtreeCursor, over a dual depth-first walk.
#   julia --project=prototype/env prototype/usecases/climate_grids/tree_parity.jl
using GridInterface
using GridInterface.Structured
import GridInterface as GrI
import ConservativeRegridding as CR
import ConservativeRegridding.Trees as Trees
import GeometryOps as GO, GeometryOpsCore as GOC, GeoInterface as GI, Extents
import SpatialTreeInterface as STI
using GeometryOps.UnitSpherical: SphericalCap, UnitSphericalPoint

# --- CR's test grids, built as test/trees/grids.jl builds them ---------------
make_lonlat_point_matrix(nx, ny) = [(lon, lat) for lon in range(-180, 180, length=nx + 1), lat in range(-90, 90, length=ny + 1)]
function make_polygon_matrix(points)
    nx, ny = size(points) .- 1
    polygons = Matrix{GI.Polygon}(undef, nx, ny)
    for i in 1:nx, j in 1:ny
        polygons[i, j] = GI.Polygon([GI.LinearRing([points[i, j], points[i+1, j], points[i+1, j+1], points[i, j+1], points[i, j]])])
    end
    return polygons
end
const SPH_LONS = collect(range(0.0, 40.0, length=5)); const SPH_LATS = collect(range(0.0, 30.0, length=4))
const SPH_PTS = GO.UnitSphereFromGeographic().([(lon, lat) for lon in SPH_LONS, lat in SPH_LATS])
const WARPED_SPH_PTS = [GO.UnitSphereFromGeographic()((-50.0 + 100.0 * u, -20.0 + 40.0 * v + 30.0 * sinpi(4.0 * u)))
                        for u in range(0.0, 1.0; length=9), v in range(0.0, 1.0; length=5)]
const GLOBAL_SPH_LONS = collect(range(-180.0, 180.0; length=9)); const GLOBAL_SPH_LATS = collect(range(-90.0, 90.0; length=5))
const GLOBAL_SPH_PTS = [GO.UnitSphereFromGeographic()((lon, lat)) for lon in GLOBAL_SPH_LONS, lat in GLOBAL_SPH_LATS]
f32(p) = UnitSphericalPoint(Float32(p[1]), Float32(p[2]), Float32(p[3]))

cases = Pair{String,Tuple{Any,Any}}[]
for (nx, ny) in ((16, 16), (13, 17), (3, 5))
    pts = make_lonlat_point_matrix(nx, ny)
    polys = make_polygon_matrix(pts)
    lons = collect(range(-180, 180, length=nx + 1)); lats = collect(range(-90, 90, length=ny + 1))
    push!(cases, "ExplicitPolygonGrid $(nx)×$(ny) (Planar)" => (Trees.ExplicitPolygonGrid(polys), Structured.ExplicitPolygonGrid(polys)))
    push!(cases, "CellBasedGrid $(nx)×$(ny) (Planar)" => (Trees.CellBasedGrid(pts), Structured.CellBasedGrid(pts)))
    push!(cases, "RegularGrid $(nx)×$(ny) (Planar)" => (Trees.RegularGrid(lons, lats), Structured.RegularGrid(lons, lats)))
end
push!(cases, "RegularGrid 4×3 regional (Spherical)" => (Trees.RegularGrid(GOC.Spherical(), SPH_LONS, SPH_LATS), Structured.RegularGrid(GOC.Spherical(), SPH_LONS, SPH_LATS)))
push!(cases, "CellBasedGrid 4×3 regional (Spherical)" => (Trees.CellBasedGrid(SPH_PTS), Structured.CellBasedGrid(SPH_PTS)))
push!(cases, "CellBasedGrid 8×4 warped (Spherical)" => (Trees.CellBasedGrid(WARPED_SPH_PTS), Structured.CellBasedGrid(WARPED_SPH_PTS)))
push!(cases, "RegularGrid 8×4 global (Spherical)" => (Trees.RegularGrid(GOC.Spherical(), GLOBAL_SPH_LONS, GLOBAL_SPH_LATS), Structured.RegularGrid(GOC.Spherical(), GLOBAL_SPH_LONS, GLOBAL_SPH_LATS)))
push!(cases, "CellBasedGrid 8×4 global (Spherical)" => (Trees.CellBasedGrid(GLOBAL_SPH_PTS), Structured.CellBasedGrid(GLOBAL_SPH_PTS)))
push!(cases, "CellBasedGrid 4×3 regional Float32 (Spherical)" => (Trees.CellBasedGrid(f32.(SPH_PTS)), Structured.CellBasedGrid(f32.(SPH_PTS))))
push!(cases, "RegularGrid 8×4 global Float32 (Spherical)" => (Trees.RegularGrid(GOC.Spherical(), Float32.(GLOBAL_SPH_LONS), Float32.(GLOBAL_SPH_LATS)),
    Structured.RegularGrid(GOC.Spherical(), Float32.(GLOBAL_SPH_LONS), Float32.(GLOBAL_SPH_LATS))))
push!(cases, "CellBasedGrid 36×18 lon/lat (Spherical)" => (Trees.CellBasedGrid(GOC.Spherical(), [GO.UnitSphereFromGeographic()((lo, la)) for lo in range(-180, 180; length=37), la in range(-90, 90; length=19)]),
    Structured.CellBasedGrid(GOC.Spherical(), [GO.UnitSphereFromGeographic()((lo, la)) for lo in range(-180, 180; length=37), la in range(-90, 90; length=19)])))
# the planar unit-square polygon grids (#65), the wide 1x1 cap cell, the authalic lon/lat grids
function make_planar_grid(nx, ny)
    polys = Matrix{GI.Polygon}(undef, nx, ny)
    for j in 1:ny, i in 1:nx
        x0, x1 = (i - 1) / nx, i / nx; y0, y1 = (j - 1) / ny, j / ny
        polys[i, j] = GI.Polygon([GI.LinearRing([(x0, y0), (x1, y0), (x1, y1), (x0, y1), (x0, y0)])])
    end
    return polys
end
for (nx, ny) in ((2, 2), (4, 3))
    polys = make_planar_grid(nx, ny)
    push!(cases, "ExplicitPolygonGrid $(nx)×$(ny) unit square (Planar)" => (Trees.ExplicitPolygonGrid(GOC.Planar(), polys), Structured.ExplicitPolygonGrid(GOC.Planar(), polys)))
end
let (sw, se, ne, nw) = GO.UnitSphereFromGeographic().(((-110.0, -70.0), (30.0, -40.0), (140.0, 20.0), (-170.0, 80.0)))
    pts = reshape(typeof(sw)[sw, se, nw, ne], 2, 2)
    push!(cases, "CellBasedGrid 1×1 wide cap (Spherical)" => (Trees.CellBasedGrid(pts), Structured.CellBasedGrid(pts)))
end
const AUTHALIC = GOC.Spherical(; radius=6371007.180918474)
authalic_corners(nlon, nlat) = [GO.UnitSphereFromGeographic()((lon, lat)) for lon in range(-180, 180; length=nlon), lat in range(-90, 90; length=nlat)]
for (nlon, nlat) in ((9, 5), (7, 4))
    pts = authalic_corners(nlon, nlat)
    push!(cases, "CellBasedGrid $(nlon-1)×$(nlat-1) authalic (Spherical)" => (Trees.CellBasedGrid(AUTHALIC, pts), Structured.CellBasedGrid(AUTHALIC, pts)))
end

# --- the dual walk ----------------------------------------------------------
exact(a::Extents.Extent, b::Extents.Extent) = a == b
exact(a::SphericalCap, b::SphericalCap) = a.point == b.point && a.radius == b.radius
exact(a, b) = false
approx(a::Extents.Extent, b::Extents.Extent) = all(isapprox.(a.X, b.X)) && all(isapprox.(a.Y, b.Y))
approx(a::SphericalCap, b::SphericalCap) = isapprox(a.point, b.point; atol=1e-12) && isapprox(a.radius, b.radius; rtol=1e-12)
approx(a, b) = false
struct Tally
    nodes::Base.RefValue{Int}; leaves::Base.RefValue{Int}
    ext_exact::Base.RefValue{Int}; ext_approx::Base.RefValue{Int}; ext_bad::Base.RefValue{Int}
    shape_bad::Base.RefValue{Int}; leaf_bad::Base.RefValue{Int}
end
Tally() = Tally(Ref(0), Ref(0), Ref(0), Ref(0), Ref(0), Ref(0), Ref(0))
function dual_dfs!(t::Tally, a, b)
    t.nodes[] += 1
    ea, eb = STI.node_extent(a), STI.node_extent(b)
    exact(ea, eb) ? (t.ext_exact[] += 1) : approx(ea, eb) ? (t.ext_approx[] += 1) : (t.ext_bad[] += 1)
    la, lb = STI.isleaf(a), STI.isleaf(b)
    if la != lb
        t.shape_bad[] += 1
        return
    end
    if la
        t.leaves[] += 1
        ia = [(i, e) for (i, e) in STI.child_indices_extents(a)]
        ib = [(i, e) for (i, e) in STI.child_indices_extents(b)]
        ok = length(ia) == length(ib) && first.(ia) == first.(ib) && all(exact(x[2], y[2]) || approx(x[2], y[2]) for (x, y) in zip(ia, ib))
        ok || (t.leaf_bad[] += 1)
        return
    end
    na, nb = STI.nchild(a), STI.nchild(b)
    if na != nb
        t.shape_bad[] += 1
        return
    end
    for k in 1:na
        dual_dfs!(t, STI.getchild(a, k), STI.getchild(b, k))
    end
end

println(rpad("grid", 50), rpad("nodes", 8), rpad("leaves", 8), rpad("extents ==", 12), rpad("extents ≈", 12), rpad("extent ≠", 10), rpad("shape ≠", 9), "leaves ≠")
allok = true
for (name, (g_cr, g_gi)) in cases
    q_cr = Trees.TopDownQuadtreeCursor(g_cr)
    q_gi = Structured.TopDownQuadtreeCursor(g_gi, 1; leafsize=(2, 2))
    t = Tally()
    dual_dfs!(t, q_cr, q_gi)
    ok = t.ext_bad[] == 0 && t.shape_bad[] == 0 && t.leaf_bad[] == 0
    global allok &= ok
    println(rpad(name, 50), rpad(t.nodes[], 8), rpad(t.leaves[], 8), rpad(t.ext_exact[], 12), rpad(t.ext_approx[], 12),
        rpad(t.ext_bad[], 10), rpad(t.shape_bad[], 9), t.leaf_bad[], ok ? "" : "   <-- differs")
end
# the same walk with the sweep leaf size, so GI's ForSweeps tree matches CR's cursor at (4, 4)
for (name, (g_cr, g_gi)) in filter(c -> occursin("36×18", c.first) || occursin("Float32 (Spherical)", c.first) && occursin("global", c.first), cases)
    t = Tally()
    dual_dfs!(t, Trees.TopDownQuadtreeCursor(g_cr; leafsize=(4, 4)), treeify(g_gi; leaf_capacity=ForSweeps()))
    ok = t.ext_bad[] == 0 && t.shape_bad[] == 0 && t.leaf_bad[] == 0
    global allok &= ok
    println(rpad(name * " leafsize (4,4)", 50), rpad(t.nodes[], 8), rpad(t.leaves[], 8), rpad(t.ext_exact[], 12), rpad(t.ext_approx[], 12),
        rpad(t.ext_bad[], 10), rpad(t.shape_bad[], 9), t.leaf_bad[], ok ? "" : "   <-- differs")
end
println("all equal: ", allok)
