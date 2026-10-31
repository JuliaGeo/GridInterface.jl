```@meta
CurrentModule = GridInterface
```

# Vertices and the dual

A vertex is a point where three or more cells meet, or a corner on the outer boundary of a regional grid. Identity has two layers: typed keys ([`VertexKey`](@ref), named locally by a grid with analytic vertices) and the dense numbering `1:nvertices(grid)`, the provider's own or derived.

```@docs
Conformity
Conforming
NonConforming
conformity
VertexKey
vertices_are_analytic
vertex_numbering_is_analytic
nvertices
cell_vertices
vertex_cells
vertex_point
vertex_isopen
vertex_key
vertex_index
cell_vertex_keys
vertexat
DualGrid
dual
dual_ring_cells
VertexTable
vertex_table
vertex_numbering
_vertex_precomputation
```
