# Tiny hand-written grids used by the unit tests and the conformance runs.
import GridInterface as GrI
using GridInterface
import GeometryOps as GO
import GeometryOpsCore as GOC
using GeometryOps.UnitSpherical: UnitSphericalPoint, SphericalCap
import Extents
import GeoInterface as GI

"A fixed-size ring with no heap storage."
struct SRing{N,P} <: AbstractVector{P}
    pts::NTuple{N,P}
end
Base.size(::SRing{N}) where {N} = (N,)
Base.IndexStyle(::Type{<:SRing}) = IndexLinear()
Base.@propagate_inbounds Base.getindex(r::SRing, i::Int) = r.pts[i]

# ---------------------------------------------------------------------------
# 1. The minimal planar grid: nx × ny unit squares, column-major. Only the
#    required set: ncells, best_manifold, cellshape, cell_boundary.
# ---------------------------------------------------------------------------
struct TinyPlanar <: GrI.AbstractGrid
    nx::Int
    ny::Int
end
GrI.ncells(g::TinyPlanar) = g.nx * g.ny
GOC.best_manifold(::TinyPlanar) = GOC.Planar()
GrI.cellshape(::TinyPlanar) = GrI.ConvexCells()
function GrI.cell_boundary(g::TinyPlanar, i::Integer)
    1 <= i <= GrI.ncells(g) || throw(BoundsError(g, i))
    x = Float64((i - 1) % g.nx); y = Float64((i - 1) ÷ g.nx)
    return SRing(((x, y), (x + 1, y), (x + 1, y + 1), (x, y + 1)))
end

# ---------------------------------------------------------------------------
# 2. A polygon bag: GeneralCells, overlapping (two overlapping squares, a
#    square with a hole, a two-part multipolygon).
# ---------------------------------------------------------------------------
struct PolygonBag <: GrI.AbstractGrid
    geoms::Vector{Any}
end
function PolygonBag()
    sq(x0, y0, s) = [(x0, y0), (x0 + s, y0), (x0 + s, y0 + s), (x0, y0 + s), (x0, y0)]
    holey = GI.Polygon([GI.LinearRing(sq(10.0, 0.0, 4.0)), GI.LinearRing(reverse(sq(11.0, 1.0, 2.0)))])
    multi = GI.MultiPolygon([GI.Polygon([GI.LinearRing(sq(20.0, 0.0, 1.0))]), GI.Polygon([GI.LinearRing(sq(22.0, 0.0, 1.0))])])
    PolygonBag(Any[GI.Polygon([GI.LinearRing(sq(0.0, 0.0, 2.0))]), GI.Polygon([GI.LinearRing(sq(1.0, 1.0, 2.0))]), holey, multi])
end
GrI.ncells(b::PolygonBag) = length(b.geoms)
GOC.best_manifold(::PolygonBag) = GOC.Planar()
GrI.tessellation(::PolygonBag) = GrI.Overlapping()
GrI.cell_geometry(b::PolygonBag, i::Integer) = b.geoms[i]

# ---------------------------------------------------------------------------
# 3. A tiny spherical grid: gnomonic cube sphere, n × n cells per face. Edges
#    are great-circle arcs, so the cells are convex. Required set only.
# ---------------------------------------------------------------------------
struct TinyCubeSphere <: GrI.AbstractGrid
    n::Int
end
const _FACES = (((1, 0, 0), (0, 1, 0), (0, 0, 1)), ((-1, 0, 0), (0, 0, 1), (0, 1, 0)),
                ((0, 1, 0), (0, 0, 1), (1, 0, 0)), ((0, -1, 0), (1, 0, 0), (0, 0, 1)),
                ((0, 0, 1), (1, 0, 0), (0, 1, 0)), ((0, 0, -1), (0, 1, 0), (1, 0, 0)))
GrI.ncells(g::TinyCubeSphere) = 6 * g.n^2
GOC.best_manifold(::TinyCubeSphere) = GOC.Spherical(; radius=1.0)
GrI.cellshape(::TinyCubeSphere) = GrI.ConvexCells()
@inline function _cubept(f, u, v)
    N, U, V = _FACES[f]
    x = N[1] + u * U[1] + v * V[1]; y = N[2] + u * U[2] + v * V[2]; z = N[3] + u * U[3] + v * V[3]
    r = sqrt(x^2 + y^2 + z^2)
    return UnitSphericalPoint(x / r, y / r, z / r)
end
function GrI.cell_boundary(g::TinyCubeSphere, i::Integer)
    1 <= i <= GrI.ncells(g) || throw(BoundsError(g, i))
    n = g.n
    f = (i - 1) ÷ n^2 + 1
    k = (i - 1) % n^2
    a = k % n; b = k ÷ n
    u0 = -1 + 2a / n; u1 = -1 + 2(a + 1) / n; v0 = -1 + 2b / n; v1 = -1 + 2(b + 1) / n
    return SRing((_cubept(f, u0, v0), _cubept(f, u1, v0), _cubept(f, u1, v1), _cubept(f, u0, v1)))
end

# ---------------------------------------------------------------------------
# 4. A planar hierarchical system: the quadtree of the unit square. Level l has
#    4^l cells in Morton order; position = code + 1. Sorted subtrees,
#    congruent refinement, analytic location and neighbours.
# ---------------------------------------------------------------------------
struct QuadCell <: GrI.AbstractCellIndex
    level::Int32
    code::Int64
end
Base.isless(a::QuadCell, b::QuadCell) = isless((a.level, a.code), (b.level, b.code))
GrI.level(c::QuadCell) = Int(c.level)

struct QuadSystem <: GrI.AbstractHierarchicalGridSystem
    maxlevel::Int
end
struct QuadLevelGrid <: GrI.AbstractGrid
    sys::QuadSystem
    level::Int
end

GrI.levels(s::QuadSystem) = 0:s.maxlevel
GrI.rootcells(::QuadSystem) = (QuadCell(0, 0),)
Base.parent(::QuadSystem, c::QuadCell) =
    c.level == 0 ? throw(ArgumentError("root has no parent")) : QuadCell(c.level - 1, c.code >> 2)
function GrI.children(s::QuadSystem, c::QuadCell)
    c.level < s.maxlevel || throw(ArgumentError("no children below maxlevel"))
    b = c.code << 2
    return (QuadCell(c.level + 1, b), QuadCell(c.level + 1, b + 1), QuadCell(c.level + 1, b + 2), QuadCell(c.level + 1, b + 3))
end
function GrI.levelgrid(s::QuadSystem, l::Integer)
    0 <= l <= s.maxlevel || throw(ArgumentError("level $l outside 0:$(s.maxlevel)"))
    return QuadLevelGrid(s, Int(l))
end
GrI.has_sorted_subtrees(::QuadSystem) = true
GrI.has_congruent_refinement(::QuadSystem) = true
function GrI.descendant_range(::QuadSystem, c::QuadCell, l::Integer)
    d = 2 * (l - c.level)
    d >= 0 || throw(ArgumentError("level $l above the cell"))
    return Int(c.code << d) + 1:Int((c.code + 1) << d)
end
GrI.node_extent(::QuadSystem, c::QuadCell) = _quadbox(c)

# Morton code <-> (ix, iy)
function _demorton(code::Int64)
    ix = 0; iy = 0
    for b in 0:30
        ix |= ((code >> (2b)) & 1) << b
        iy |= ((code >> (2b + 1)) & 1) << b
    end
    return ix, iy
end
function _morton(ix::Int, iy::Int)
    code = Int64(0)
    for b in 0:30
        code |= Int64((ix >> b) & 1) << (2b)
        code |= Int64((iy >> b) & 1) << (2b + 1)
    end
    return code
end
function _quadbox(c::QuadCell)
    ix, iy = _demorton(c.code); h = 1.0 / (1 << c.level)
    return Extents.Extent(X=(ix * h, (ix + 1) * h), Y=(iy * h, (iy + 1) * h))
end

GrI.system(g::QuadLevelGrid) = g.sys
GrI.level(g::QuadLevelGrid) = g.level
GrI.ncells(g::QuadLevelGrid) = 1 << (2 * g.level)
GOC.best_manifold(::QuadLevelGrid) = GOC.Planar()
GrI.cellshape(::QuadLevelGrid) = GrI.ConvexCells()
GrI.extent_is_expensive(::QuadLevelGrid) = false
function GrI.cellindex(g::QuadLevelGrid, i::Integer)
    1 <= i <= GrI.ncells(g) || throw(BoundsError(g, i))
    return QuadCell(g.level, i - 1)
end
GrI.localindex(g::QuadLevelGrid, c::QuadCell) =
    (c.level == g.level && 0 <= c.code < GrI.ncells(g)) ? Int(c.code) + 1 : nothing
function GrI.cell_boundary(g::QuadLevelGrid, i::Integer)
    e = GrI.cell_extent(g, i)
    return SRing(((e.X[1], e.Y[1]), (e.X[2], e.Y[1]), (e.X[2], e.Y[2]), (e.X[1], e.Y[2])))
end
GrI.cell_extent(g::QuadLevelGrid, i::Integer) = _quadbox(GrI.cellindex(g, i))
GrI.cell_centroid(g::QuadLevelGrid, i::Integer) = (e = GrI.cell_extent(g, i); ((e.X[1] + e.X[2]) / 2, (e.Y[1] + e.Y[2]) / 2))
GrI.cell_area(g::QuadLevelGrid, i::Integer) = 1.0 / GrI.ncells(g)
GrI.location_is_analytic(::QuadLevelGrid) = true
function GrI.cellat(g::QuadLevelGrid, p)
    x, y = Float64(GI.x(p)), Float64(GI.y(p))
    (0 <= x <= 1 && 0 <= y <= 1) || return nothing
    k = 1 << g.level
    ix = min(floor(Int, x * k), k - 1); iy = min(floor(Int, y * k), k - 1)
    return Int(_morton(ix, iy)) + 1
end
GrI.neighbors_are_analytic(::QuadLevelGrid) = true
GrI.maxneighbors(::QuadLevelGrid, ::GrI.Vertex) = 8
GrI.maxneighbors(::QuadLevelGrid, ::GrI.Edge) = 4
GrI.winding(::QuadLevelGrid, ::GrI.Connectivity) = GrI.CounterClockwise()
# counter-clockwise from the neighbour across edge 1 (the bottom edge), aligned with the corners
const _VOFF = ((0, -1), (1, -1), (1, 0), (1, 1), (0, 1), (-1, 1), (-1, 0), (-1, -1))
const _EOFF = ((0, -1), (1, 0), (0, 1), (-1, 0))
function GrI.neighbors(g::QuadLevelGrid, i::Integer; connectivity::GrI.Connectivity=GrI.Vertex())
    ix, iy = _demorton(Int64(i - 1)); k = 1 << g.level
    offs = connectivity isa GrI.Edge ? _EOFF : _VOFF
    out = Int[]
    for (dx, dy) in offs
        jx, jy = ix + dx, iy + dy
        (0 <= jx < k && 0 <= jy < k) && push!(out, Int(_morton(jx, jy)) + 1)
    end
    return out
end
GrI.chunks(g::QuadLevelGrid) = g.level >= 1 ? GrI.hierarchy_chunks(g, 1) : nothing

# ---------------------------------------------------------------------------
# 5. The four-method spherical grid: lon/lat cells with great-circle edges.
#    Pole cells repeat the pole point (their corners collapse to three).
# ---------------------------------------------------------------------------
struct MinLonLat <: GrI.AbstractGrid
    nlon::Int
    nlat::Int
end
GrI.ncells(g::MinLonLat) = g.nlon * g.nlat
GOC.best_manifold(::MinLonLat) = GOC.Spherical(; radius=1.0)
GrI.cellshape(::MinLonLat) = GrI.SimpleCells()
function GrI.cell_boundary(g::MinLonLat, i::Integer)
    1 <= i <= GrI.ncells(g) || throw(BoundsError(g, i))
    a = (i - 1) % g.nlon; b = (i - 1) ÷ g.nlon
    x0 = -180 + 360a / g.nlon; x1 = -180 + 360(a + 1) / g.nlon
    y0 = -90 + 180b / g.nlat; y1 = -90 + 180(b + 1) / g.nlat
    f = GO.UnitSphereFromGeographic()
    return SRing((f((x0, y0)), f((x1, y0)), f((x1, y1)), f((x0, y1))))
end

# The same grid handing out lon/lat tuples on a spherical manifold: refused.
struct TupleLonLat <: GrI.AbstractGrid end
GrI.ncells(::TupleLonLat) = 4
GOC.best_manifold(::TupleLonLat) = GOC.Spherical(; radius=1.0)
GrI.cellshape(::TupleLonLat) = GrI.SimpleCells()
GrI.cell_boundary(::TupleLonLat, i::Integer) =
    SRing(((10.0i, 0.0), (10.0i + 10, 0.0), (10.0i + 10, 10.0), (10.0i, 10.0)))
