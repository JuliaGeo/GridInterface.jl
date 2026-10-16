# ===========================================================================
# Transparent wrappers: a grid plus something the caller holds for it (a tree,
# a chunking, a vertex table). Positions are unchanged, so every verb forwards
# to the parent once and a wrapper overrides the one or two verbs it changes.
# `Subset`, `DualGrid` and `MultiLevelSet` remap positions: not wrappers.
#
# Caches are owned by the caller through these: GridInterface holds no global
# state, and a fallback that needs a tree or a vertex table builds it per call
# unless the grid it is given is the wrapper that stores one.
# ===========================================================================

"""
    abstract type AbstractWrappedGrid

Supertype of GridInterface's transparent wrappers ([`GridAndTree`](@ref),
[`ChunkedGrid`](@ref), [`GridAndVertices`](@ref)): a grid with the same
positions as `parent(g)` plus something the caller holds for it. Every verb in
[`WRAPPED_VERBS`](@ref) forwards to `parent(g)`; a subtype defines
`Base.parent` and overrides what it changes. Providers may subtype it for
their own pass-through types but never have to.
"""
abstract type AbstractWrappedGrid end
isgrid(::AbstractWrappedGrid) = true
Base.broadcastable(g::AbstractWrappedGrid) = Ref(g)

"""
    WRAPPED_VERBS

The verbs an [`AbstractWrappedGrid`](@ref) forwards, as `(verb, argument
types after the grid)`; keyword arguments pass through. One entry per method
shape of the verb, so that each forward is at least as specific as the generic
method it shadows (a single varargs forward would be ambiguous with them).
The conformance suite checks a wrapper against this list.
"""
const WRAPPED_VERBS = (
    # identity
    (:ncells, ()), (:cellindextypes, ()), (:cellindex, (Integer,)), (:cellindex, (Integer, Type)),
    (:localindex, (Integer,)), (:localindex, (Any,)), (:globalindex, (Any,)),
    (:level, ()), (:system, ()), (:gridrelation, ()), (:parentpositions, ()),
    # geometry and location
    (:best_manifold, ()), (:cellshape, ()), (:edgekind, ()), (:tessellation, ()), (:isglobal, ()),
    (:cellsize, ()), (:cell_geometry_is_expensive, ()), (:extent_is_expensive, ()),
    (:location_is_analytic, ()),
    (:cell_boundary, (Integer,)), (:cell_geometry, (Integer,)), (:cell_corners, (Integer,)),
    (:cell_centroid, (Integer,)), (:cell_area, (Integer,)), (:cell_extent, (Integer,)),
    (:cell_contains, (Integer, Any)), (:cellat, (Any,)), (:cellsat, (Any,)),
    (:query, (QueryPredicate,)), (:query, (Any,)),
    # tree
    (:leaf_capacity, (LeafHint,)), (:leaf_capacity, (Balanced,)), (:leaf_capacity, (ForSweeps,)),
    (:leaf_capacity, (ForQueries,)), (:treeify, (LeafHint,)), (:treeify, (Integer,)),
    (:subtree, (AbstractVector{<:Integer},)), (:subtree, (Subset,)), (:grid_extent, ()),
    # subsets and chunks
    (:chunks, ()), (:chunks, (Integer,)), (:chunkat, (Integer,)),
    # neighbours
    (:neighbors_are_analytic, ()), (:maxneighbors, ()), (:maxneighbors, (Connectivity,)),
    (:maxneighbors, (Integer,)), (:maxneighbors, (Integer, Connectivity)),
    (:maxring, (Integer,)), (:maxring, (Integer, Connectivity)),
    (:winding, ()), (:winding, (Connectivity,)),
    (:neighbors, (Integer,)), (:neighbors, (Integer, Integer)),
    (:neighbors, (Integer, Index{Local})), (:neighbors, (Integer, Index{Global})),
    (:neighbors, (Integer, Index{Type{T}} where {T})),
    (:ring, (Integer, Integer)),
    (:halo, ()), (:halo, (Integer,)), (:halo, (AbstractVector{<:Integer},)),
    (:halo, (AbstractVector{<:Integer}, Integer)),
    (:border, ()), (:border, (AbstractVector{<:Integer},)),
    (:interior, ()), (:interior, (AbstractVector{<:Integer},)),
    # dual
    (:conformity, ()), (:vertices_are_analytic, ()), (:vertex_numbering_is_analytic, ()),
    (:nvertices, ()), (:dual, ()),
    (:cell_vertices, (Integer,)), (:cell_vertex_keys, (Integer,)),
    (:vertex_cells, (Integer,)), (:vertex_cells, (VertexKey,)),
    (:vertex_point, (Integer,)), (:vertex_point, (VertexKey,)),
    (:vertex_isopen, (Integer,)), (:vertex_isopen, (VertexKey,)),
    (:vertex_key, (Integer,)), (:vertex_key, (Integer, Integer)), (:vertex_index, (VertexKey,)),
    (:vertex_table, ()), (:vertex_numbering, ()), (:vertexat, (Any,)),
)

_fwd_args(Ts) = [Symbol(:a, k) for k in eachindex(Ts)]
_fwd_sig(Ts) = [:($a::$T) for (a, T) in zip(_fwd_args(Ts), Ts)]

for (f, Ts) in WRAPPED_VERBS
    args = _fwd_args(Ts)
    @eval $f(g::AbstractWrappedGrid, $(_fwd_sig(Ts)...); kw...) = $f(parent(g), $(args...); kw...)
end
# the two verbs whose grid is not the first argument
reindex(::Type{T}, g::AbstractWrappedGrid, c) where {T} = reindex(T, parent(g), c)
cell_boundary!(buf::AbstractVector, g::AbstractWrappedGrid, i::Integer) = cell_boundary!(buf, parent(g), i)

# ---------------------------------------------------------------------------
# A grid paired with a tree: the unit engines consume.
# ---------------------------------------------------------------------------

"""
    GridAndTree(grid, tree)
    GridAndTree(grid; leaf_capacity = Balanced())   # = GridAndTree(grid, treeify(grid; leaf_capacity))

What flows into an engine (a conservative regridder, a chunked plan): a grid
paired with a tree. Engines ask the tree traversal questions and the grid
everything else (cells, polygons, areas, manifold), so the tree itself needs
no reference back to its grid. The pair forwards every grid verb to `grid`,
returns `tree` from `treeify` whatever the hint or capacity, and forwards
`split_weight` and `subtree_indices` to `tree`. Point location and queries on
the pair descend the stored tree.

Wrap once to reuse an expensive tree (an R-tree, a cap-cached cursor, a
provider's own tree) wherever a grid is accepted: this is the only tree cache
in GridInterface. With `tree = subtree(grid, cells)` the leaves emit `grid`'s
positions.

A tree does not survive chunking: `chunks` and `subtree` of the pair answer for
the inner grid. Cache per chunk instead, as `GridAndTree(chunk, treeify(chunk))`.
"""
struct GridAndTree{G,T} <: AbstractWrappedGrid
    grid::G
    tree::T
end
GridAndTree(grid; leaf_capacity=Balanced()) = GridAndTree(grid, treeify(grid; leaf_capacity))

Base.parent(g::GridAndTree) = g.grid
treeify(g::GridAndTree; kw...) = g.tree
treeify(g::GridAndTree, ::LeafHint) = g.tree
treeify(g::GridAndTree, ::Integer) = g.tree
split_weight(g::GridAndTree) = split_weight(g.tree)
subtree_indices(g::GridAndTree) = subtree_indices(g.tree)
grid_extent(g::GridAndTree) = STI.node_extent(g.tree)

# The grid's own analytic answer if it has one; otherwise the tree-backed
# fallback on the pair itself, whose `treeify` is the stored tree.
cellat(g::GridAndTree, p) = location_is_analytic(g.grid) ? cellat(g.grid, p) :
    fallback_cellat(best_manifold(g), g, check_point(best_manifold(g), p))
cellsat(g::GridAndTree, p) = fallback_cellsat(best_manifold(g), g, check_point(best_manifold(g), p))
query(g::GridAndTree, pred::QueryPredicate) = fallback_query(best_manifold(g), g, pred)
query(g::GridAndTree, region) = query(g, Intersects(region))

# ---------------------------------------------------------------------------
# A grid with a chunking other than its storage partition.
# ---------------------------------------------------------------------------

"""
    ChunkedGrid(grid, chunks)

`grid` with `chunks` answering `chunks` (and [`chunkat`](@ref) consistently);
everything else forwards. What [`mockchunks`](@ref) returns.
"""
struct ChunkedGrid{G,C} <: AbstractWrappedGrid
    grid::G
    chunks::C
end
Base.parent(g::ChunkedGrid) = g.grid
chunks(g::ChunkedGrid) = g.chunks
chunkat(g::ChunkedGrid, i::Integer) = _chunkat(g.chunks, g, Int(i))

"""
    mockchunks(grid, chunks) -> ChunkedGrid
    mockchunks(grid, LevelChunks(l))

A grid whose `chunks` is the given partition (a vector of subsets of `grid`,
or `nothing`) instead of its storage partition, for a consumer that wants to
work in different blocks than the data is stored in. [`LevelChunks`](@ref)`(l)`
names the level-`l` blocks of a hierarchical grid. Mirrors DiskArrays' verb of
the same name; composes with subsets and [`GridAndTree`](@ref) in any order.
"""
mockchunks(grid, c) = ChunkedGrid(grid, c)
mockchunks(grid, c::LevelChunks) = ChunkedGrid(grid, chunks(grid, c.level))

# ---------------------------------------------------------------------------
# A grid paired with its vertex precomputation: the dual's primal.
# ---------------------------------------------------------------------------

"""
    GridAndVertices(grid)
    GridAndVertices(grid, table)

`grid` paired with the precomputation its dense vertex numbering needs: the
geometric [`VertexTable`](@ref) (the fallback's merged corners), the derived
numbering of a grid with analytic vertices, or `nothing` for a grid whose
numbering is analytic ([`vertex_numbering_is_analytic`](@ref)). The dense
vertex verbs (`nvertices`, `cell_vertices`, `vertex_cells`, `vertex_point`,
`vertex_isopen`, `vertex_key`, `vertex_index`, `cell_vertex_keys`, `vertexat`)
read the stored table; everything else forwards. Build once where those verbs
are called repeatedly: on a plain grid each call rebuilds the table, which is
O(ncells). [`dual`](@ref) holds the same table.
"""
struct GridAndVertices{G,T} <: AbstractWrappedGrid
    grid::G
    table::T       # VertexTable (geometric vertices), VertexNumbering (analytic vertices), or nothing
end
GridAndVertices(grid) = GridAndVertices(grid, _vertex_precomputation(grid))

Base.parent(g::GridAndVertices) = g.grid
vertex_table(g::GridAndVertices{G,<:VertexTable}) where {G} = g.table
vertex_numbering(g::GridAndVertices{G,VertexNumbering}) where {G} = g.table
dual(g::GridAndVertices) = DualGrid(g.grid, g.table)

# The dense verbs run their generic bodies on the pair, which read the stored
# table through the two methods above, instead of forwarding. A grid whose
# numbering is analytic answers them itself (its methods, not the generic
# bodies, which an `invoke` on the pair would bypass): forward to it.
for (f, Ts) in ((:nvertices, ()), (:cell_vertices, (Integer,)), (:cell_vertex_keys, (Integer,)),
                (:vertex_cells, (Integer,)), (:vertex_point, (Integer,)),
                (:vertex_isopen, (Integer,)), (:vertex_isopen, (VertexKey,)),
                (:vertex_key, (Integer,)), (:vertex_index, (VertexKey,)), (:vertexat, (Any,)))
    args = _fwd_args(Ts)
    @eval $f(g::GridAndVertices, $(_fwd_sig(Ts)...)) = vertex_numbering_is_analytic(g.grid) ?
        $f(g.grid, $(args...)) : invoke($f, Tuple{Any,$(Ts...)}, g, $(args...))
end
# The keyed forms are the provider's on a grid with analytic vertices.
vertex_key(g::GridAndVertices, i::Integer, k::Integer) = vertices_are_analytic(g) ?
    vertex_key(g.grid, i, k) : invoke(vertex_key, Tuple{Any,Integer,Integer}, g, i, k)
vertex_cells(g::GridAndVertices, key::VertexKey) = vertices_are_analytic(g) ?
    vertex_cells(g.grid, key) : invoke(vertex_cells, Tuple{Any,VertexKey}, g, key)
