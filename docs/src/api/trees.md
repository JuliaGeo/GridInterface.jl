```@meta
CurrentModule = GridInterface
```

# Trees

The contracts are described in [Trees and restriction](../trees.md).

```@docs
LeafHint
Balanced
ForSweeps
ForQueries
treeify
leaf_capacity
SubtreeCapacity
own_tree
subtree
grid_extent
split_weight
subtree_indices
position_span
span_range
span_mayhold
span_inside
LatticeBlock
MaskedCursor
```

## The hierarchy cursor

The generic tree of a level grid, a subset of one, or a multi-level set, written against the hierarchy primitives.

```@docs
HierarchyCursor
SelectionCursor
RootIndex
root_index
hierarchy_tree
LevelFrontier
SubsetFrontier
GlobalKeys
LeafEntries
```

## Queries

```@docs
QueryPredicate
Intersects
Within
Contains
CentroidIn
query
query_candidates
query_region
```
