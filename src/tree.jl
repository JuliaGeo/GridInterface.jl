# ===========================================================================
# Tree family: you get a SpatialTreeInterface tree from a grid.
#
# Consumers ask the tree only traversal questions (children, extents, the
# positions its leaves emit). Cell count, polygons and manifold are asked of the
# grid. What flows into an engine is therefore a grid paired with a tree:
# `GridAndTree(grid, treeify(grid))` (see wrappers.jl).
# ===========================================================================

"""
    abstract type LeafHint

What a tree is for; [`treeify`](@ref) asks the grid for the leaf size that
suits it through [`leaf_capacity`](@ref)`(grid, hint)`:

| hint | leaf size (core default) | for |
|---|---|---|
| [`Balanced`](@ref)`()` | 4 | mixed use (the default) |
| [`ForSweeps`](@ref)`()` | 8 | dual-tree sweeps (a regridder): bucketed leaves |
| [`ForQueries`](@ref)`()` | 1 | point location and queries: the finest tree in O(1) |
"""
abstract type LeafHint end
"The default intent: a middle ground between sweeps and queries."
struct Balanced <: LeafHint end
"Dual-tree sweeps: bucketed leaves (one-cell leaves made CR's sweep visit 27x more node pairs)."
struct ForSweeps <: LeafHint end
"Point location and queries: the finest tree the provider has in O(1)."
struct ForQueries <: LeafHint end

"""
    treeify(grid; leaf_capacity = Balanced())
    treeify(grid, hint::LeafHint)
    treeify(grid, n::Integer)

A SpatialTreeInterface tree whose leaves emit `grid`'s positions, each exactly
once. The keyword names the intent (a [`LeafHint`](@ref)) or a leaf size; the
core turns an intent into a number with [`leaf_capacity`](@ref)`(grid, hint)`
and calls the integer form, which is the provider's seam: a provider defines
`treeify(g, n::Integer)` (its one tree builder) and `leaf_capacity(g, hint)`
for the intents it tunes; a provider with trees of different kinds may catch
an intent directly with `treeify(g, ::ForQueries)`. A
[`GridAndTree`](@ref) ignores both and returns its stored tree.

Generic integer form:

- a level grid of a hierarchical system, or a subset of one (`system(grid)` and
  `level(grid)` answer), gets the generic hierarchy cursor: O(1) to build when
  the system has sorted subtrees, with `n` members per leaf and the system's
  [`root_index`](@ref) above its roots;
- anything else gets an STR R-tree over [`cell_extent`](@ref) with node
  capacity `n` (extension), O(n log n) to build: wrap it once in
  [`GridAndTree`](@ref) to reuse it.
"""
treeify(grid; leaf_capacity=Balanced()) = treeify(grid, leaf_capacity)
treeify(grid, hint::LeafHint) = treeify(grid, leaf_capacity(grid, hint))
treeify(grid, n::Integer) = _default_treeify(grid, n)

"""
    SubtreeCapacity(n) <: Integer

The leaf size [`subtree`](@ref) hands `treeify(grid, n)` when it asks for the
grid's own tree to restrict. A provider's `treeify(g, n::Integer)` reads it as
the integer `n` (it converts, promotes and compares as one); the generic
`treeify` answers `nothing` where it would otherwise build the R-tree over the
whole grid, since a subtree then builds its R-tree over just the cells instead.
So the two generic paths share one builder (`_default_treeify`), and the
exception costs no dispatch of its own.
"""
struct SubtreeCapacity <: Integer
    n::Int
end
(::Type{T})(c::SubtreeCapacity) where {T<:Integer} = T(c.n)
Base.promote_rule(::Type{SubtreeCapacity}, ::Type{T}) where {T<:Integer} = promote_type(Int, T)

# The one generic tree builder, shared by `treeify` and `subtree`: a level
# grid's hierarchy cursor, else the R-tree over the whole grid. Only the R-tree
# costs O(ncells), and a `SubtreeCapacity` turns it into `nothing`.
function _default_treeify(grid, n::Integer)
    if system(grid) !== nothing && level(grid) !== nothing
        return hierarchy_tree(grid; bucket=Int(n))
    end
    return _whole_rtree(grid, n)
end
_whole_rtree(grid, n::Integer) = fallback_treeify(best_manifold(grid), grid, Int(n))
_whole_rtree(grid, ::SubtreeCapacity) = nothing

"""
    own_tree(grid, n::Integer) -> tree or nothing

The tree `subtree` restricts: `treeify(grid, n)`, unless that would be
GridInterface's R-tree over the whole grid (`nothing`; see
[`SubtreeCapacity`](@ref)).
"""
own_tree(grid, n::Integer) = treeify(grid, SubtreeCapacity(Int(n)))

"""
    leaf_capacity(grid, hint::LeafHint) -> Int

How many cells a leaf of `treeify(grid, hint)` holds: the generic hierarchy
cursor's bucket, the R-tree's node capacity, a provider's own bucket. Core
defaults: `Balanced` 4, `ForSweeps` 8, `ForQueries` 1; a provider overrides
them per intent. Never read inside a node: only `treeify` consumes it.
"""
leaf_capacity(grid, ::Balanced) = 4
leaf_capacity(grid, ::ForSweeps) = 8
leaf_capacity(grid, ::ForQueries) = 1

"""
    subtree(grid, cells; leaf_capacity = Balanced())
    subtree(grid, s::Subset; leaf_capacity = Balanced())

A tree over the cells at positions `cells` of `grid` whose leaves emit
**`grid`'s** positions (parent positions): the grid's own tree,
`treeify(grid; leaf_capacity)`, restricted to `cells` by a
[`MaskedCursor`](@ref), whatever package built that tree. The restriction
keeps `cells` lazy (a range stays a range), prunes children by
[`position_span`](@ref), filters leaves by membership and starts at the
deepest node that alone may hold members, so a chunk roots at its own node in
O(depth). `cells` must be strictly ascending (an unsorted vector is sorted
first); a [`Subset`](@ref) of this very grid restricts to its positions.

Dispatch on the collection lets a typed one build in O(1): a raster tile type
starts the provider's quadtree at the tile, a multi-level set's chunk its
cursor at its root. The exception is a grid with no tree of its own, whose
`treeify` would be GridInterface's R-tree over the whole grid (see
[`SubtreeCapacity`](@ref)): its subtree is the R-tree over just `cells`.

The keyword is [`treeify`](@ref)'s: a [`LeafHint`](@ref), turned into a leaf
size by [`leaf_capacity`](@ref)`(grid, hint)`, or the leaf size itself.
"""
function subtree(grid, cells::AbstractVector{<:Integer}; leaf_capacity=Balanced())
    n = _capacity(grid, leaf_capacity)
    t = own_tree(grid, n)
    t === nothing && return fallback_subtree(best_manifold(grid), grid, cells, n)
    return MaskedCursor{false}(t, _ascending(cells))
end

# `treeify`'s reading of the keyword, for the `subtree` methods
_capacity(grid, hint::LeafHint) = leaf_capacity(grid, hint)
_capacity(grid, n::Integer) = Int(n)

"""
    split_weight(node) -> Real

SpatialTreeInterface's weight of a node (its number of entries, or an
estimate that never grows going down the tree), by which a dual traversal
chooses the side to open. GridInterface adds methods for its own trees and
wrappers: a hierarchy cursor's window size, a quadtree block's cell count, a
[`GridAndTree`](@ref)'s stored tree's weight, a [`MaskedCursor`](@ref)'s
span overlap with its members.
"""
split_weight

"""
    subtree_indices(node) -> AbstractVector{<:Integer}

SpatialTreeInterface's leaf indices under `node`, in traversal order, without
computing extents. GridInterface's trees answer in O(1) where the leaves
occupy a contiguous window (a hierarchy cursor's keys, a quadtree block); a
[`MaskedCursor`](@ref) answers the members under the node.
"""
subtree_indices

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
on the sphere). Candidates come from descending
`treeify(grid; leaf_capacity = ForQueries())` (a [`GridAndTree`](@ref)'s stored
tree) with the region's covering extent; the exact test runs on each candidate
(extension). A provider with a dedicated query engine overrides it. Nothing is
cached: a caller who queries repeatedly wraps the grid once in `GridAndTree`.
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
    STI.depth_first_search(e -> Extents.intersects(e, q), treeify(grid; leaf_capacity=ForQueries())) do i
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
