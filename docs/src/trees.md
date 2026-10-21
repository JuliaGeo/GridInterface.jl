```@meta
CurrentModule = GridInterface
```

# Trees and restriction

A grid's tree is a [SpatialTreeInterface](https://github.com/JuliaGeo/GeometryOps.jl/tree/as/split-unitsphericalgeometrycore/SpatialTreeInterface) tree whose leaves emit the grid's positions, each exactly once. Consumers ask the tree traversal questions only (children, node extents, a leaf's `(position, extent)` entries) and the grid everything else, so what flows into an engine is a grid paired with a tree: [`GridAndTree`](@ref)`(grid, treeify(grid))`.

## `treeify`: the grid's own tree

[`treeify`](@ref) turns an intent or a leaf size into one call of the provider's seam:

| call | what it does |
|---|---|
| `treeify(grid)`, `treeify(grid, hint::LeafHint)` | asks [`leaf_capacity`](@ref)`(grid, hint)` for the leaf size, then calls the integer form |
| `treeify(grid, n::Integer)` | the provider's one tree builder |

The generic integer form builds a [`HierarchyCursor`](@ref) for a level grid of a hierarchical system (O(1) to build when the system has sorted subtrees), and an STR R-tree over [`cell_extent`](@ref) for anything else (O(n log n); wrap it in [`GridAndTree`](@ref) to reuse it).

A provider types its leaf size as `n::Integer`. [`subtree`](@ref) passes a [`SubtreeCapacity`](@ref) through `treeify`, and a method typed `n::Int` would be skipped for it; the conformance law "subtree restricts treeify(grid, n)'s own tree" enforces this.

## `subtree`: the grid's own tree over some of its cells

[`subtree`](@ref)`(grid, cells; leaf_capacity)` restricts the grid's own tree to the cells at positions `cells`, whatever package built that tree, with leaves emitting `grid`'s positions. It is one generic function over [`MaskedCursor`](@ref); providers implement `treeify` and, on their node types, [`position_span`](@ref), and never override `subtree`. The restriction

1. takes the grid's own tree, [`own_tree`](@ref)`(grid, n)`;
2. descends to the tight root, the deepest node that alone may hold members, so a chunk roots at its own node in O(depth);
3. drops children whose `position_span` misses every member;
4. filters leaf entries to members.

`position_span(node)` is a bound on the parent positions beneath a node:

| answer | meaning |
|---|---|
| `UnitRange{Int}` | the positions lie in this run |
| [`LatticeBlock`](@ref) | `rows` runs of `width` positions in a column-major lattice, exact against a range of members |
| a provider's own type | any value answering [`span_range`](@ref), [`span_mayhold`](@ref) and [`span_inside`](@ref) (DiscreteGlobalGrids' `TileRows`) |
| `nothing` | the tree cannot say cheaply (an R-tree); the child is kept |

A bound, never the exact set: a loose bound costs visits, never correctness. Every tree GridInterface builds answers.

### `SubtreeCapacity`

[`SubtreeCapacity`](@ref)` <: Integer` is the leaf size `subtree` hands `treeify(grid, n)`. A provider's `treeify(g, n::Integer)` reads it as the integer `n`. The generic `treeify` answers `nothing` with it where it would otherwise build the R-tree over the whole grid; `subtree` then builds its R-tree over just `cells` (`fallback_subtree`). This is the one case where `subtree` restricts no tree: the grid has none of its own.

### A subset's own tree

`treeify(s::Subset, n)` is the same restriction with leaves emitting the subset's own positions (`1:ncells(s)`), the `L = true` form of `MaskedCursor`; `subtree(grid, s)` for a `Subset` of that grid emits the parent's. A parent with no tree of its own leaves the subset its R-tree.

## Node extent cost

`SpatialTreeInterface.node_extent_is_expensive(::Type{Node})` declares, per node type, whether `node_extent` computes (a union of caps, a perimeter walk) or reads. A dual-tree engine reads it to decide which side's extents to hold. GridInterface's cursors answer `true`, the planar `RegularGrid` quadtree `false`, and `MaskedCursor` forwards to the node type it wraps. The cost of an extent is a property of the tree, so there is no grid-level trait for it ([`cell_geometry_is_expensive`](@ref) is the grid's trait for polygons).

## Leaf entries

`child_indices_extents(leaf)` is built once per leaf visit. A bucketed leaf (a leaf size above one, the sweep shape) collects its entries, so a dual traversal that reads the inner side's leaf once per opposing cell reads stored extents rather than recomputing them `n` times. A one-cell leaf (the query shape) stays lazy ([`LeafEntries`](@ref)). The structured quadtree builds its entries in an `NTuple` with no heap storage ([`Structured.LeafEntries`](@ref)).

## Queries

[`query`](@ref)`(grid, predicate)` descends `treeify(grid; leaf_capacity = ForQueries())` by extent and runs the exact predicate on each candidate ([`Intersects`](@ref), [`Within`](@ref), [`Contains`](@ref), [`CentroidIn`](@ref)). Nothing is cached: a caller who queries repeatedly wraps the grid once in `GridAndTree`.
