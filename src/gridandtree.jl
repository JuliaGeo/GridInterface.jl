# ===========================================================================
# A grid paired with a tree: the unit engines consume.
# ===========================================================================

"""
    GridAndTree(grid, tree)
    GridAndTree(grid)                       # = GridAndTree(grid, treeify(grid))

What flows into an engine (a conservative regridder, a chunked plan): a grid
paired with a tree. Engines ask the tree traversal questions and the grid
everything else (cells, polygons, areas, manifold), so the tree itself needs
no reference back to its grid. The pair forwards every grid verb to `grid`,
returns `tree` from `treeify`, and forwards [`split_weight`](@ref) to `tree`.

Wrap once to reuse an expensive tree (an R-tree, a cap-cached cursor, a
provider's own tree) wherever a grid is accepted. With
`tree = subtree(grid, cells)` the leaves emit `grid`'s positions.

A tree does not survive chunking: `chunks` and `subtree` of the pair answer for
the inner grid. Cache per chunk instead, as `GridAndTree(chunk, treeify(chunk))`.
"""
struct GridAndTree{G,T} <: AbstractGrid
    grid::G
    tree::T
end
GridAndTree(grid) = GridAndTree(grid, treeify(grid))

treeify(g::GridAndTree) = g.tree
Base.parent(g::GridAndTree) = g.grid
split_weight(g::GridAndTree) = split_weight(g.tree)

# Every grid verb forwards. Kept as explicit lists so the conformance suite can
# check them against the verbs it knows.
const FORWARDED_TRAITS = (
    :ncells, :cellindextypes, :system, :level, :gridrelation, :cellshape, :tessellation,
    :conformity, :location_is_analytic, :neighbors_are_analytic, :extent_is_expensive,
    :nvertices, :chunks, :maxneighbors, :winding, :edgekind, :cell_geometry_is_expensive,
    :cellsize, :isglobal, :leaf_capacity, :vertices_are_analytic,
)
for f in FORWARDED_TRAITS
    @eval $f(g::GridAndTree) = $f(g.grid)
end
GOC.best_manifold(g::GridAndTree) = best_manifold(g.grid)
grid_extent(g::GridAndTree) = STI.node_extent(g.tree)
maxneighbors(g::GridAndTree, c::Connectivity) = maxneighbors(g.grid, c)
maxneighbors(g::GridAndTree, k::Integer, c::Connectivity=Vertex()) = maxneighbors(g.grid, k, c)
maxring(g::GridAndTree, k::Integer, c::Connectivity=Vertex()) = maxring(g.grid, k, c)
winding(g::GridAndTree, c::Connectivity) = winding(g.grid, c)
chunks(g::GridAndTree, l::Integer) = chunks(g.grid, l)
chunkat(g::GridAndTree, i::Integer) = chunkat(g.grid, i)

const FORWARDED_CELL_VERBS = (:cell_boundary, :cell_geometry, :cell_centroid, :cell_area, :cell_extent,
    :cell_corners, :cell_vertices, :cell_vertex_keys)
for f in FORWARDED_CELL_VERBS
    @eval $f(g::GridAndTree, i::Integer) = $f(g.grid, i)
end
cell_boundary!(buf::AbstractVector, g::GridAndTree, i::Integer) = cell_boundary!(buf, g.grid, i)
cellindex(g::GridAndTree, i::Integer) = cellindex(g.grid, i)
cellindex(g::GridAndTree, i::Integer, ::Type{T}) where {T} = cellindex(g.grid, i, T)
reindex(::Type{T}, g::GridAndTree, c) where {T} = reindex(T, g.grid, c)
localindex(g::GridAndTree, c) = localindex(g.grid, c)
localindex(g::GridAndTree, c::Integer) = localindex(g.grid, c)
globalindex(g::GridAndTree, c) = globalindex(g.grid, c)
cell_contains(g::GridAndTree, i::Integer, p) = cell_contains(g.grid, i, p)
vertex_cells(g::GridAndTree, v::Integer) = vertex_cells(g.grid, v)
vertex_cells(g::GridAndTree, k::VertexKey) = vertex_cells(g.grid, k)
vertex_point(g::GridAndTree, v::Integer) = vertex_point(g.grid, v)
vertex_point(g::GridAndTree, k::VertexKey) = vertex_point(g.grid, k)
vertex_isopen(g::GridAndTree, v::Integer) = vertex_isopen(g.grid, v)
vertex_key(g::GridAndTree, i::Integer, k::Integer) = vertex_key(g.grid, i, k)
vertex_key(g::GridAndTree, v::Integer) = vertex_key(g.grid, v)
vertex_index(g::GridAndTree, k::VertexKey) = vertex_index(g.grid, k)
neighbors(g::GridAndTree, i::Integer; connectivity::Connectivity=Vertex()) = neighbors(g.grid, i; connectivity)
neighbors(g::GridAndTree, i::Integer, k::Integer; connectivity::Connectivity=Vertex()) = neighbors(g.grid, i, k; connectivity)
neighbors(g::GridAndTree, i::Integer, sel::Index{Local}; connectivity::Connectivity=Vertex()) = neighbors(g.grid, i, sel; connectivity)
neighbors(g::GridAndTree, i::Integer, sel::Index{Global}; connectivity::Connectivity=Vertex()) = neighbors(g.grid, i, sel; connectivity)
neighbors(g::GridAndTree, i::Integer, sel::Index{Type{T}}; connectivity::Connectivity=Vertex()) where {T} =
    neighbors(g.grid, i, sel; connectivity)
ring(g::GridAndTree, i::Integer, k::Integer; connectivity::Connectivity=Vertex()) = ring(g.grid, i, k; connectivity)
halo(g::GridAndTree, positions::AbstractVector{<:Integer}, k::Integer=1; connectivity::Connectivity=Vertex()) =
    halo(g.grid, positions, k; connectivity)
border(g::GridAndTree, positions::AbstractVector{<:Integer}; connectivity::Connectivity=Vertex()) =
    border(g.grid, positions; connectivity)
interior(g::GridAndTree, positions::AbstractVector{<:Integer}; connectivity::Connectivity=Vertex()) =
    interior(g.grid, positions; connectivity)
subtree(g::GridAndTree, cells::AbstractVector{<:Integer}) = subtree(g.grid, cells)
subtree(g::GridAndTree, s::Subset) = subtree(g.grid, s)
dual(g::GridAndTree) = dual(g.grid)

# Point location and queries: the grid's own analytic answer if it has one;
# otherwise the tree-backed fallback over the *stored* tree.
cellat(g::GridAndTree, p) = location_is_analytic(g.grid) ? cellat(g.grid, p) :
    fallback_cellat(best_manifold(g), g, check_point(best_manifold(g), p))
cellsat(g::GridAndTree, p) = fallback_cellsat(best_manifold(g), g, check_point(best_manifold(g), p))
