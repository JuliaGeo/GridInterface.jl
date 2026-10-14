# ===========================================================================
# Cell geometry, the shape trait, the tessellation law, point location.
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
"""
abstract type CellShape end
abstract type Polygonal <: CellShape end
abstract type SingleRing <: Polygonal end
struct ConvexCells <: SingleRing end
struct SimpleCells <: SingleRing end
struct GeneralCells <: Polygonal end

"""
    cellshape(grid) -> CellShape

Declared per grid, not per cell. Default `GeneralCells()`. For a hierarchical
system the class is declared once and holds for every level grid.
"""
cellshape(grid) = GeneralCells()

"""
    best_manifold(grid)

The manifold the grid lives on (`Planar()`, `Spherical()`, ...). **Required.**
This is `GeometryOpsCore.best_manifold`, the same binding ConservativeRegridding
calls, so a grid that defines it is understood by the regridder too.
"""
best_manifold

# ---------------------------------------------------------------------------
# A zero-copy closed ring and a one-element vector, so the generic
# `cell_geometry` of a single-ring grid allocates nothing beyond what the
# grid's own `cell_boundary` allocates.
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

# ---------------------------------------------------------------------------
# The geometry verbs. Position-first; typed ids are accepted and resolved
# through `localindex` (which the PRD expects to be analytic).
# Providers should type the position argument `::Integer` (or `::Int`).
# ---------------------------------------------------------------------------

@noinline _refuse_boundary(grid) = throw(ArgumentError(
    "cell_boundary exists only for SingleRing grids; $(typeof(grid)) declares $(cellshape(grid))"))
@noinline _missing_boundary(grid) = throw(ArgumentError(
    "$(typeof(grid)) declares $(cellshape(grid)) but defines no `cell_boundary`"))
@noinline _missing_geometry(grid) = throw(ArgumentError(
    "$(typeof(grid)) declares GeneralCells and must define `cell_geometry`"))

"""
    cell_boundary(grid, i) -> AbstractVector of points

The exterior ring of cell `i`: counter-clockwise, **unclosed**. Exists only for
`SingleRing` grids (refused otherwise). A single-ring grid may implement this
alone and get [`cell_geometry`](@ref) for free.
"""
cell_boundary(grid, i::Integer) = _cell_boundary(cellshape(grid), grid, i)
_cell_boundary(::SingleRing, grid, i) = _missing_boundary(grid)
_cell_boundary(::CellShape, grid, i) = _refuse_boundary(grid)

"""
    cell_geometry(grid, i)

Cell `i` as a GeoInterface polygonal geometry (polygon, possibly with holes, or
multipolygon) whose vertices lie on the manifold and whose edges are its
geodesics. The primitive of the geometry family. Generic method for
`SingleRing` grids: the closed ring of [`cell_boundary`](@ref).
"""
cell_geometry(grid, i::Integer) = _cell_geometry(cellshape(grid), grid, i)
_cell_geometry(::SingleRing, grid, i) = ring_polygon(cell_boundary(grid, i))
_cell_geometry(::CellShape, grid, i) = _missing_geometry(grid)

"""
    cell_centroid(grid, i)

A point strictly inside cell `i`. Fallback (extension): from the geometry.
"""
cell_centroid(grid, i::Integer) = fallback_centroid(best_manifold(grid), grid, i)

"""
    cell_area(grid, i)

Area of cell `i` on the grid's manifold (steradians on the unit sphere).
Fallback (extension): `GeometryOps.area` of the geometry.
"""
cell_area(grid, i::Integer) = fallback_area(best_manifold(grid), grid, i)

"""
    cell_extent(grid, i)

The per-cell extent the tree family prunes with: a box (`Extents.Extent`) on
the plane, a spherical cap on the sphere (by convention). It must contain the
cell's geometry. Fallback (extension): from the geometry.
"""
cell_extent(grid, i::Integer) = fallback_extent(best_manifold(grid), grid, i)

# Typed-id forms: one `localindex`, then the position form.
for f in (:cell_boundary, :cell_geometry, :cell_centroid, :cell_area, :cell_extent)
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

Whether a point lies in the interior of at most one cell. `Tessellating()`
(the default and the law) or `Overlapping()` (a polygon bag such as country
outlines, or a multi-level set on a system whose children do not tile their
parent). On an overlapping grid [`cellat`](@ref) returns the lowest matching
position unless the grid answers analytically.
"""
abstract type Tessellation end
struct Tessellating <: Tessellation end
struct Overlapping <: Tessellation end

tessellation(grid) = Tessellating()

# ---------------------------------------------------------------------------
# Point location.
# ---------------------------------------------------------------------------

"""
    location_is_analytic(grid) -> Bool

Whether [`cellat`](@ref) answers without a tree descent (a DGGS inverse
projection, a raster's index arithmetic). Consumers that walk from a seed cell
check it and use the tree otherwise. Default `false`.
"""
location_is_analytic(grid) = false

"""
    cellat(grid, p) -> Union{Int,Nothing}

Position of the cell containing point `p` (a unit-sphere point on spherical
grids), or `nothing`. Always available: an analytic grid overrides it; any other
grid gets the tree-backed fallback (descend `treeify(grid)` by extent, then an
interior test), which returns the **lowest** containing position. That is the
tie rule on shared boundaries and on overlapping grids.
"""
cellat(grid, p) = fallback_cellat(best_manifold(grid), grid, p)

"""
    cell_contains(grid, i, p) -> Bool

Closed interior test of point `p` against cell `i` (boundary counts as inside).
Fallback (extension): half-plane tests on `ConvexCells`, winding number
otherwise.
"""
cell_contains(grid, i::Integer, p) = fallback_contains(best_manifold(grid), grid, i, p)
