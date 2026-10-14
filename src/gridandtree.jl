# ===========================================================================
# Caching a tree is a transparent wrapper, not a keyword.
# ===========================================================================

"""
    GridAndTree(grid, tree)

A grid that forwards every grid verb to `grid` and returns `tree` from
`treeify`. Wrap once to reuse an expensive tree (an R-tree, a cap-cached
cursor) wherever a grid is accepted. With `tree = subtree(grid, cells)` the
leaves emit `grid`'s positions, which is how a consumer regrids a chunk in the
global index space.
"""
struct GridAndTree{G,T} <: AbstractGrid
    grid::G
    tree::T
end

treeify(g::GridAndTree) = g.tree
Base.parent(g::GridAndTree) = g.grid

# Every grid verb forwards. Kept as one explicit list so the conformance suite
# can check it against the verbs it knows.
const FORWARDED_VERBS = (
    :ncells, :cellindextypes, :system, :level, :gridrelation, :cellshape, :tessellation,
    :conformity, :location_is_analytic, :neighbors_are_analytic, :extent_is_expensive,
    :nvertices, :chunks, :maxneighbors, :winding,
)
for f in FORWARDED_VERBS
    @eval $f(g::GridAndTree) = $f(g.grid)
end
GOC.best_manifold(g::GridAndTree) = best_manifold(g.grid)
grid_extent(g::GridAndTree) = STI.node_extent(g.tree)
maxneighbors(g::GridAndTree, c::Connectivity) = maxneighbors(g.grid, c)
winding(g::GridAndTree, c::Connectivity) = winding(g.grid, c)

for f in (:cell_boundary, :cell_geometry, :cell_centroid, :cell_area, :cell_extent,
          :cell_vertices)
    @eval $f(g::GridAndTree, i::Integer) = $f(g.grid, i)
end
cellindex(g::GridAndTree, i::Integer) = cellindex(g.grid, i)
cellindex(g::GridAndTree, i::Integer, ::Type{T}) where {T} = cellindex(g.grid, i, T)
reindex(::Type{T}, g::GridAndTree, c) where {T} = reindex(T, g.grid, c)
localindex(g::GridAndTree, c) = localindex(g.grid, c)
localindex(g::GridAndTree, c::Integer) = localindex(g.grid, c)
globalindex(g::GridAndTree, c) = globalindex(g.grid, c)
cell_contains(g::GridAndTree, i::Integer, p) = cell_contains(g.grid, i, p)
vertex_cells(g::GridAndTree, v::Integer) = vertex_cells(g.grid, v)
vertex_point(g::GridAndTree, v::Integer) = vertex_point(g.grid, v)
neighbors(g::GridAndTree, i::Integer; connectivity::Connectivity=Vertex()) = neighbors(g.grid, i; connectivity)
ring(g::GridAndTree, i::Integer, k::Integer; connectivity::Connectivity=Vertex()) = ring(g.grid, i, k; connectivity)
subtree(g::GridAndTree, cells::AbstractVector{<:Integer}) = subtree(g.grid, cells)
dual(g::GridAndTree) = dual(g.grid)

# Point location: the grid's own analytic answer if it has one; otherwise the
# tree-backed fallback over the *stored* tree.
cellat(g::GridAndTree, p) = location_is_analytic(g.grid) ? cellat(g.grid, p) :
    fallback_cellat(best_manifold(g), g, p)
