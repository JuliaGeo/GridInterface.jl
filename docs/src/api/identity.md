```@meta
CurrentModule = GridInterface
```

# Identity

The grid owns a dense order, positions `1:ncells(grid)`, which every verb takes and returns. Typed ids ([`AbstractCellIndex`](@ref)) are a layer on top, reached through [`cellindex`](@ref) and [`localindex`](@ref).

```@docs
AbstractGrid
isgrid
AbstractCellIndex
iscellindex
ncells
cellindex
cellindextypes
reindex
localindex
globalindex
level
Index
Local
Global
OUTSIDE
toindex
GridRelation
NoParent
ParentGrid
ParentSystem
gridrelation
system
```
