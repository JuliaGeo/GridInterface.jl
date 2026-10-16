# ===========================================================================
# Tree family: you get a SpatialTreeInterface tree from a grid.
#
# Consumers ask the tree only traversal questions (children, extents, the
# positions its leaves emit). Cell count, polygons and manifold are asked of the
# grid. What flows into an engine is therefore a grid paired with a tree:
# `GridAndTree(grid, treeify(grid))` (see gridandtree.jl).
# ===========================================================================

"""
    treeify(grid)

A SpatialTreeInterface tree whose leaves emit `grid`'s positions, each exactly
once. A provider's own `treeify` always wins (DGG's cursors, CR's quadtree, a
raster's index arithmetic). Generic method:

- a level grid of a hierarchical system, or a subset of one (`system(grid)` and
  `level(grid)` answer), gets the generic hierarchy cursor: O(1) to build when
  the system has sorted subtrees, with [`leaf_capacity`](@ref) members per leaf
  and the system's [`root_index`](@ref) above its roots;
- anything else gets an STR R-tree over [`cell_extent`](@ref) (extension).
"""
treeify(grid) = _default_treeify(grid)

function _default_treeify(grid)
    if system(grid) !== nothing && level(grid) !== nothing
        return hierarchy_tree(grid; bucket=_bucket(grid))
    end
    return fallback_treeify(best_manifold(grid), grid)
end

"""
    leaf_capacity(grid) -> Union{Int,Nothing}

Hint: how many cells a leaf of the grid's tree should hold (the generic
hierarchy cursor's bucket; the R-tree's node capacity). `nothing` (default)
lets the tree choose (1 for the hierarchy cursor, 16 for the R-tree).
"""
leaf_capacity(grid) = nothing
_bucket(grid) = (b = leaf_capacity(grid); b === nothing ? 1 : Int(b))

"""
    subtree(grid, cells)

A tree over the cells at positions `cells` of `grid` whose leaves emit
**`grid`'s** positions (parent positions). Dispatches on the type of `cells`, so
a typed collection builds in O(1): a [`Subset`](@ref) chunk with a root starts
the hierarchy cursor at that root; a raster tile type starts the provider's
quadtree at the tile. Plain vectors (sorted, unique) get the hierarchy cursor
restricted to them on level grids, an R-tree over those cells otherwise.
"""
subtree(grid, cells::AbstractVector{<:Integer}) = _default_subtree(grid, cells)

function _default_subtree(grid, cells)
    if system(grid) !== nothing && level(grid) !== nothing && gridrelation(grid) isa NoParent
        return hierarchy_subtree(grid, cells; bucket=_bucket(grid))
    end
    return fallback_subtree(best_manifold(grid), grid, cells)
end

"""
    extent_is_expensive(grid) -> Bool

Traversal hint: whether [`cell_extent`](@ref) computes rather than reads.
Default `true`. (The node-level counterpart is `STI.node_extent_is_expensive`;
the polygon counterpart is [`cell_geometry_is_expensive`](@ref).)
"""
extent_is_expensive(grid) = true

# `split_weight` and `subtree_indices` are SpatialTreeInterface's; GridInterface
# adds methods for its own trees and wrappers.

"""
    grid_extent(grid)

An extent covering every cell of `grid` (an `Extents.Extent` box on the plane, a
`SphericalCap` on the sphere). Cheap by contract for chunks (stored or O(1)).
Generic: the root extent of `treeify(grid)`. Named `grid_extent` rather than
`extent` to avoid colliding with `Extents.extent`.
"""
grid_extent(grid) = STI.node_extent(treeify(grid))

# ---------------------------------------------------------------------------
# Queries: candidates by extent, then an exact predicate.
# ---------------------------------------------------------------------------

"""
    abstract type QueryPredicate

What [`query`](@ref) tests each candidate cell against. The region `r` is a
GeoInterface polygon or multipolygon (with `UnitSphericalPoint` vertices on
the sphere), an `Extents.Extent`, or (sphere) a `SphericalCap`:

| predicate | a cell is selected when |
|---|---|
| [`Intersects`](@ref)`(r)` | the closed cell and `r` share at least a point |
| [`Within`](@ref)`(r)` | the cell lies inside `r` (boundaries may touch) |
| [`Contains`](@ref)`(r)` | the cell contains all of `r` |
| [`CentroidIn`](@ref)`(r)` | [`cell_centroid`](@ref) lies in closed `r` |

**A lon/lat `Extents.Extent{(:X,:Y)}` on the sphere** is read as DGG reads it:
the box boundary sampled every 2° or finer along its parallels and meridians,
joined by great-circle arcs, as a polygon. Candidates come from the cap that
covers it ([`query_region`](@ref)); the predicate is then exact against that
polygon. On the plane an `Extent` is the rectangle.
"""
abstract type QueryPredicate end
"Select cells sharing at least a point with the region (closed sets)."
struct Intersects{G} <: QueryPredicate
    region::G
end
"Select cells lying inside the region (touching its boundary allowed)."
struct Within{G} <: QueryPredicate
    region::G
end
"Select cells that contain the whole region."
struct Contains{G} <: QueryPredicate
    region::G
end
"Select cells whose centroid lies in the closed region."
struct CentroidIn{G} <: QueryPredicate
    region::G
end

"""
    query(grid, pred::QueryPredicate) -> Vector{Int}
    query(grid, region)                    # = query(grid, Intersects(region))

Ascending positions of the cells satisfying `pred` (see
[`QueryPredicate`](@ref) for the predicates and the reading of a lon/lat extent
on the sphere). Candidates come from descending the cached `treeify(grid)`
with the region's covering extent; the exact test runs on each candidate
(extension). A provider with a dedicated query engine overrides it.
"""
query(grid, pred::QueryPredicate) = fallback_query(best_manifold(grid), grid, pred)
query(grid, region) = query(grid, Intersects(region))

"""
    query_candidates(grid, region) -> Vector{Int}

Ascending positions of the cells whose [`cell_extent`](@ref) intersects the
covering extent of `region` (no exact test: a superset of every predicate's
answer).
"""
function query_candidates(grid, region)
    m = best_manifold(grid)
    q = query_region(m, region)
    out = Int[]
    STI.depth_first_search(e -> Extents.intersects(e, q), treeify(grid)) do i
        push!(out, i)
        nothing
    end
    return sort!(out)
end

"""
    query_region(manifold, region)

The covering extent of `region` in the grid's own extent type. Identity for an
extent of the right type; the extension maps a lon/lat `Extents.Extent` to a
covering `SphericalCap` on the sphere and a polygon to its covering extent.
"""
query_region(m, region) = fallback_region_extent(m, region)
query_region(::Any, ext::Extents.Extent) = ext
query_region(::Any, ext::SphericalCap) = ext
