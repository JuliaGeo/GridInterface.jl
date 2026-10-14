# ===========================================================================
# Cell geometry, the shape and edge traits, the tessellation law, point location.
# ===========================================================================

"""
    abstract type CellShape

What a grid promises about every one of its cells; a hierarchy so consumers
dispatch on "at least this strong".

    CellShape
    └─ Polygonal                 polygon with holes, or multipolygon
       ├─ GeneralCells           anything polygonal (the default)
       └─ SingleRing             one exterior ring, no holes, non-self-intersecting
          ├─ SimpleCells         ring may be concave
          └─ ConvexCells         ring is convex and inside an open hemisphere

A consumer refuses a shape weaker than its engine needs, up front: the
spherical conservative engine needs `ConvexCells` on the destination (clip)
side only; general clipping on the sphere is not available, so a planar view of
a spherical grid is the pattern for zonal statistics over country polygons.
"""
abstract type CellShape end
"Abstract: polygonal cells (with holes, or multipolygons, allowed)."
abstract type Polygonal <: CellShape end
"Abstract: every cell is one exterior ring with no holes; [`cell_boundary`](@ref) exists."
abstract type SingleRing <: Polygonal end
"Every cell is a convex single ring inside an open hemisphere."
struct ConvexCells <: SingleRing end
"Every cell is a simple (possibly concave) single ring."
struct SimpleCells <: SingleRing end
"Cells may have holes or several parts; only [`cell_geometry`](@ref) exists."
struct GeneralCells <: Polygonal end

"""
    cellshape(grid) -> CellShape

**Required** (one of the four required methods; see [`GridInterface`](@ref)).
Declared per grid, not per cell. The default is `GeneralCells()`, which refuses
`cell_boundary`, so a single-ring grid must declare `SimpleCells()` or
`ConvexCells()`. For a hierarchical system the class is declared once and holds
for every level grid.
"""
cellshape(grid) = GeneralCells()

"""
    abstract type EdgeKind

How the points of [`cell_boundary`](@ref) relate to the cell's true edges,
answered by [`edgekind`](@ref):

- [`GeodesicEdges`](@ref): consecutive boundary points are joined by geodesics
  (great circles on the sphere, straight lines on the plane) and the boundary
  *is* the cell. Boundary points are the corners (up to repeats).
- [`CurvedEdges`](@ref): the true edges are curves (HEALPix, ISEA, A5 chart
  edges); the boundary is a densified approximation and
  [`cell_corners`](@ref) are the true corners. A plotter trusts the points; a
  consumer that needs exact edges asks the provider.
"""
abstract type EdgeKind end
"Boundary points joined by geodesics are the exact cell."
struct GeodesicEdges <: EdgeKind end
"True edges are curves; the boundary is densified."
struct CurvedEdges <: EdgeKind end

"""
    edgekind(grid) -> EdgeKind

Default `GeodesicEdges()`.
"""
edgekind(grid) = GeodesicEdges()

"""
    best_manifold(grid)

**Required.** The manifold the grid lives on (`Planar()`, `Spherical(radius = r)`,
...). This is `GeometryOpsCore.best_manifold`, the binding ConservativeRegridding
calls.

On a spherical (or ellipsoidal) manifold every point is a `UnitSphericalPoint`
and every extent a `SphericalCap`; the core refuses other types. The radius is
the grid's **area unit**: [`cell_area`](@ref) is in units of `radius²`, and
the cells of a global grid sum to `4π radius²`. Consumers taking two grids
never reconcile radii implicitly: two grids with different radii need an
explicit manifold (ConservativeRegridding's three-argument `Regridder`), and
bit-parity is defined per manifold.
"""
best_manifold

# ---------------------------------------------------------------------------
# Zero-copy ring views.
# ---------------------------------------------------------------------------

"""
    ClosedRing(points)

Lazy view of an unclosed ring with its first point repeated at the end.
"""
struct ClosedRing{T,V<:AbstractVector{T}} <: AbstractVector{T}
    points::V
end
Base.size(r::ClosedRing) = (length(r.points) + 1,)
Base.IndexStyle(::Type{<:ClosedRing}) = IndexLinear()
Base.@propagate_inbounds Base.getindex(r::ClosedRing, i::Int) =
    i == length(r.points) + 1 ? r.points[1] : r.points[i]

struct OneVector{T} <: AbstractVector{T}
    x::T
end
Base.size(::OneVector) = (1,)
Base.IndexStyle(::Type{<:OneVector}) = IndexLinear()
Base.@propagate_inbounds function Base.getindex(v::OneVector, i::Int)
    @boundscheck i == 1 || throw(BoundsError(v, i))
    return v.x
end

"""
    ring_polygon(boundary) -> GI.Polygon

A GeoInterface polygon over the closed view of an unclosed ring.
"""
function ring_polygon(b::AbstractVector)
    ring = GI.LinearRing(ClosedRing(b))
    return GI.Polygon(OneVector(ring))
end
ring_polygon(b) = ring_polygon(collect(b))

# coordinates of a point as a tuple of Float64 (2 on the plane, 3 on the sphere)
@inline _coords(p::UnitSphericalPoint) = (Float64(p[1]), Float64(p[2]), Float64(p[3]))
@inline _coords(p::NTuple{2,Real}) = (Float64(p[1]), Float64(p[2]))
@inline _coords(p::NTuple{3,Real}) = (Float64(p[1]), Float64(p[2]), Float64(p[3]))
@inline _coords(p) = GI.is3d(p) ? (Float64(GI.x(p)), Float64(GI.y(p)), Float64(GI.z(p))) :
    (Float64(GI.x(p)), Float64(GI.y(p)), 0.0)
@inline function _near(a, b)
    x = _coords(a); y = _coords(b)
    for k in eachindex(x)
        abs(x[k] - y[k]) <= 1e-12 * (1 + abs(x[k])) || return false
    end
    return true
end

"""
    CollapsedRing(ring)

Lazy view of an unclosed ring with repeated points collapsed: a point equal (to
1e-12 relative) to its predecessor is skipped, and so is a last point equal to
the first. The default [`cell_corners`](@ref). Length is counted once at
construction; iteration is one pass with no allocation.
"""
struct CollapsedRing{T,V<:AbstractVector{T}} <: AbstractVector{T}
    ring::V
    n::Int
end
function CollapsedRing(r::V) where {T,V<:AbstractVector{T}}
    m = length(r)
    m == 0 && return CollapsedRing{T,V}(r, 0)
    n = 1
    @inbounds for k in 2:m
        _near(r[k], r[k-1]) || (n += 1)
    end
    # a last run equal to the first point closes the ring onto itself
    k = m
    @inbounds while k > 1 && _near(r[k], r[k-1])
        k -= 1
    end
    @inbounds (k > 1 && n > 1 && _near(r[k], r[1])) && (n -= 1)
    return CollapsedRing{T,V}(r, n)
end
Base.size(c::CollapsedRing) = (c.n,)
Base.IndexStyle(::Type{<:CollapsedRing}) = IndexLinear()
function Base.iterate(c::CollapsedRing, st=(1, 0))
    k, emitted = st
    emitted >= c.n && return nothing
    r = c.ring
    @inbounds while k > 1 && k <= length(r) && _near(r[k], r[k-1])
        k += 1
    end
    return (@inbounds(r[k]), (k + 1, emitted + 1))
end
function Base.getindex(c::CollapsedRing, i::Int)
    @boundscheck 1 <= i <= c.n || throw(BoundsError(c, i))
    for (k, p) in enumerate(c)
        k == i && return p
    end
    error("unreachable")
end

# ---------------------------------------------------------------------------
# The geometry verbs. Position-first; typed ids are accepted and resolved
# through `localindex` (which the PRD expects to be analytic).
# Providers should type the position argument `::Integer` (or `::Int`).
# ---------------------------------------------------------------------------

@noinline _refuse_boundary(grid) = throw(ArgumentError(
    "cell_boundary exists only for SingleRing grids; $(typeof(grid)) declares $(cellshape(grid)). " *
    "A single-ring grid must declare `GridInterface.cellshape(::MyGrid) = SimpleCells()` (or ConvexCells())."))
@noinline _missing_boundary(grid) = throw(ArgumentError(
    "$(typeof(grid)) declares $(cellshape(grid)) but defines no `cell_boundary`"))
@noinline _missing_geometry(grid) = throw(ArgumentError(
    "$(typeof(grid)) declares GeneralCells (the default `cellshape`) and defines no `cell_geometry`. " *
    "A single-ring grid declares `cellshape(::MyGrid) = SimpleCells()` and defines `cell_boundary`."))

"""
    cell_boundary(grid, i) -> AbstractVector of points

**Required for single-ring grids** (the fourth required method, or
[`cell_geometry`](@ref)). The exterior ring of cell `i`: counter-clockwise,
**unclosed**, points on the manifold (`UnitSphericalPoint`s on the sphere).
Exists only for `SingleRing` grids (refused otherwise). Return a static ring
(`SVector`, a tuple-backed vector) to make sweeps allocation free, or define
[`cell_boundary!`](@ref).
"""
cell_boundary(grid, i::Integer) = _cell_boundary(cellshape(grid), grid, i)
_cell_boundary(::SingleRing, grid, i) = _missing_boundary(grid)
_cell_boundary(::CellShape, grid, i) = _refuse_boundary(grid)

"""
    cell_boundary!(buf::AbstractVector, grid, i) -> buf

Non-allocating ring access: `buf` resized to and filled with
[`cell_boundary`](@ref)`(grid, i)`. The generic method copies (so it allocates
whatever `cell_boundary` allocates); a provider whose rings are computed (a
densified chart boundary) overrides it to write into `buf` directly.
"""
function cell_boundary!(buf::AbstractVector, grid, i::Integer)
    b = cell_boundary(grid, i)
    resize!(buf, length(b))
    copyto!(buf, b)
    return buf
end

"""
    cell_corners(grid, i) -> AbstractVector of points

The distinct true corners of cell `i`, counter-clockwise, unclosed. Default:
[`cell_boundary`](@ref) with repeated points collapsed (a lon/lat pole cell
lists the pole once), lazily. A grid with [`CurvedEdges`](@ref) overrides it
to return the corners without the densification points. Vertex identity and
the dual are built on corners.
"""
cell_corners(grid, i::Integer) = CollapsedRing(cell_boundary(grid, i))

"""
    cell_geometry(grid, i)

**Required for general grids** (instead of [`cell_boundary`](@ref)). Cell `i`
as a GeoInterface polygonal geometry (polygon, possibly with holes, or
multipolygon) whose vertices lie on the manifold. Generic method for
`SingleRing` grids: the closed ring of [`cell_boundary`](@ref), no copy.
"""
cell_geometry(grid, i::Integer) = _cell_geometry(cellshape(grid), grid, i)
_cell_geometry(::SingleRing, grid, i) = ring_polygon(cell_boundary(grid, i))
_cell_geometry(::CellShape, grid, i) = _missing_geometry(grid)

"""
    cell_geometry_is_expensive(grid) -> Bool

Hint: whether [`cell_geometry`](@ref)/[`cell_boundary`](@ref) computes (a
densified chart boundary, a projection) rather than reads. A consumer that
visits a cell's polygon many times (per chunk pair) memoises it when this is
`true`. Default `false`. (The extent counterpart is [`extent_is_expensive`](@ref).)
"""
cell_geometry_is_expensive(grid) = false

"""
    cell_centroid(grid, i)

A point strictly inside cell `i`. Fallback (extension): from the geometry.
"""
cell_centroid(grid, i::Integer) = fallback_centroid(best_manifold(grid), grid, i)

"""
    cell_area(grid, i)

Area of cell `i` in units of the manifold's `radius²` (steradians on the unit
sphere). Fallback (extension): `GeometryOps.area` of the geometry.
"""
cell_area(grid, i::Integer) = fallback_area(best_manifold(grid), grid, i)

"""
    cell_extent(grid, i)

The per-cell extent the tree family prunes with: an `Extents.Extent` box on the
plane, a `SphericalCap` on the sphere (enforced). It must contain the cell's
geometry. Fallback (extension): from the geometry.
"""
cell_extent(grid, i::Integer) = fallback_extent(best_manifold(grid), grid, i)

"""
    cellsize(grid) -> Float64

A typical cell size, for level of detail and tolerances: the square root of the
mean cell area over up to 64 evenly spaced cells, in the manifold's length
unit (radians × radius on the sphere). Providers with an analytic answer
override it.
"""
function cellsize(grid)
    n = ncells(grid)
    n == 0 && return 0.0
    m = min(n, 64)
    s = 0.0
    for k in 1:m
        s += Float64(cell_area(grid, 1 + ((k - 1) * (n - 1)) ÷ max(m - 1, 1)))
    end
    return sqrt(s / m)
end

"""
    isglobal(grid) -> Bool

Whether the grid's cells cover the whole (spherical) manifold. A global
tessellating grid's areas sum to the manifold's area (a conformance law), and
the dual of a global conforming grid tiles the sphere. Default `false`.
"""
isglobal(grid) = false

# Typed-id forms: one `localindex`, then the position form.
for f in (:cell_boundary, :cell_geometry, :cell_centroid, :cell_area, :cell_extent, :cell_corners)
    @eval $f(grid, c::AbstractCellIndex) = $f(grid, _pos(grid, c))
end
@inline function _pos(grid, c)
    i = localindex(grid, c)
    i === nothing && throw(ArgumentError("cell $c is not a member of $(typeof(grid))"))
    return i
end

# ---------------------------------------------------------------------------
# The tessellation law.
# ---------------------------------------------------------------------------

"""
    abstract type Tessellation

Whether a point lies in the interior of at most one cell, answered by
[`tessellation`](@ref).
"""
abstract type Tessellation end
"Every point lies in the interior of at most one cell (the default, and the law)."
struct Tessellating <: Tessellation end
"Cells may overlap (a polygon bag; a multi-level set of a non-congruent system)."
struct Overlapping <: Tessellation end

"""
    tessellation(grid) -> Tessellation

`Tessellating()` (default) or `Overlapping()`. On an overlapping grid
[`cellat`](@ref) returns the lowest containing position unless the grid answers
analytically; [`cellsat`](@ref) returns every containing position.
"""
tessellation(grid) = Tessellating()

# ---------------------------------------------------------------------------
# Point location.
# ---------------------------------------------------------------------------

"""
    location_is_analytic(grid) -> Bool

Whether [`cellat`](@ref) answers without a tree descent (a DGGS inverse
projection, a raster's index arithmetic). Default `false`. An analytic answer
must agree with [`cell_contains`](@ref): `cellat(g, p) == i ⇒ cell_contains(g, i, p)`
on points near edges too (a conformance law), so a grid whose true edges are
curved densifies its boundary.
"""
location_is_analytic(grid) = false

"""
    cellat(grid, p) -> Union{Int,Nothing}

Position of a cell containing point `p` (a `UnitSphericalPoint` on spherical
grids; anything else is refused), or `nothing`. Always available: an analytic
grid overrides it; any other grid gets the tree-backed fallback (descend the
cached `treeify(grid)` by extent, then [`cell_contains`](@ref)), which returns
the **lowest** containing position: the tie rule on shared boundaries and on
overlapping grids. An analytic answer on a tessellating grid may break a
boundary tie differently; on an overlapping grid it may be any containing cell.
"""
cellat(grid, p) = (m = best_manifold(grid); fallback_cellat(m, grid, check_point(m, p)))

"""
    cellsat(grid, p) -> Vector{Int}

Ascending positions of **every** cell containing `p` (closed cells). On an
overlapping grid this is what a `sum`/`count` rasterisation or a zonal over
disputed areas needs; on a tessellating grid it is one cell, or the cells
sharing a boundary `p` lies on. Fallback (extension): tree descent plus
[`cell_contains`](@ref).
"""
cellsat(grid, p) = (m = best_manifold(grid); fallback_cellsat(m, grid, check_point(m, p)))

"""
    cell_contains(grid, i, p) -> Bool

Closed interior test of point `p` against cell `i` (boundary counts as inside).
Fallback (extension): half-plane tests on `ConvexCells`, winding number
otherwise.
"""
cell_contains(grid, i::Integer, p) =
    (m = best_manifold(grid); fallback_contains(m, grid, i, check_point(m, p)))
