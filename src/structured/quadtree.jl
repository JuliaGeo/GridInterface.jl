# ===========================================================================
# ConservativeRegridding's top-down quadtree over one face of any structured
# grid: a node is an index block (irange × jrange) of face `face`, split in
# halves per axis until both sides fit `leafsize`; a leaf emits its cells
# (i fastest) with their extents. The root tree over faces has one such
# cursor per face.
# ===========================================================================

"""
    TopDownQuadtreeCursor(grid, face = 1; leafsize = (2, 2))
    TopDownQuadtreeCursor(grid, face, irange, jrange; leafsize = (2, 2))

The quadtree over the cell block `irange × jrange` (the whole face by
default) of face `face`: CR's `TopDownQuadtreeCursor` on the face hooks. A
node's extent is [`blockextent`](@ref) of its block (by default the bounding
box or cap from the block's perimeter vertices); leaves hold at most `leafsize[1] × leafsize[2]` cells and emit
`grid`'s positions. `split_weight` is the block's cell count,
`subtree_indices` its positions in leaf order and [`block_indices`](@ref)
the same positions as a lazy ascending block.
"""
struct TopDownQuadtreeCursor{G,L,C}      # L: emit tile-local positions; C: leaf capacity
    grid::G
    face::Int
    irange::UnitRange{Int}
    jrange::UnitRange{Int}
    leafsize::NTuple{2,Int}
    i0::Int                                        # tile origin and width for local emission
    j0::Int
    ni::Int
end
function TopDownQuadtreeCursor(grid, face::Integer, irange::UnitRange{Int}, jrange::UnitRange{Int};
        leafsize::NTuple{2,Int}=(2, 2), local_positions::Bool=false)
    (leafsize[1] > 0 && leafsize[2] > 0) ||
        throw(ArgumentError("leafsize must contain positive integers; got $leafsize"))
    # one static capacity for every leaf size up to 4x4, so the intents share a type
    C = max(16, leafsize[1] * leafsize[2])
    return TopDownQuadtreeCursor{typeof(grid),local_positions,C}(grid, Int(face), irange, jrange, leafsize,
        first(irange), first(jrange), length(irange))
end
function TopDownQuadtreeCursor(grid, face::Integer=1; leafsize::NTuple{2,Int}=(2, 2))
    ni, nj = facesize(grid, face)
    return TopDownQuadtreeCursor(grid, face, 1:ni, 1:nj; leafsize)
end
_child(q::TopDownQuadtreeCursor{G,L,C}, ir, jr) where {G,L,C} =
    TopDownQuadtreeCursor{G,L,C}(q.grid, q.face, ir, jr, q.leafsize, q.i0, q.j0, q.ni)

function Base.show(io::IO, q::TopDownQuadtreeCursor)
    print(io, "TopDownQuadtreeCursor(face ", q.face, ", (", q.irange, "), (", q.jrange, "))")
end
Base.show(io::IO, ::MIME"text/plain", q::TopDownQuadtreeCursor) = show(io, q)

STI.isspatialtree(::Type{<:TopDownQuadtreeCursor}) = true
# A node's extent is `blockextent`, derived from the block's perimeter unless
# the grid has a closed form: a planar `RegularGrid`'s is four axis reads.
STI.node_extent_is_expensive(::Type{<:TopDownQuadtreeCursor}) = true
STI.node_extent_is_expensive(::Type{<:TopDownQuadtreeCursor{<:RegularGrid{<:Planar}}}) = false

STI.isleaf(q::TopDownQuadtreeCursor) =
    length(q.irange) <= q.leafsize[1] && length(q.jrange) <= q.leafsize[2]

function STI.nchild(q::TopDownQuadtreeCursor)
    i_is_leaf = length(q.irange) <= q.leafsize[1]
    j_is_leaf = length(q.jrange) <= q.leafsize[2]
    (i_is_leaf && j_is_leaf) && error("a leaf has no children")
    return (i_is_leaf || j_is_leaf) ? 2 : 4
end

@inline _lo(r::UnitRange{Int}, s) = r[1:s]
@inline _hi(r::UnitRange{Int}, s) = r[s+1:end]

function STI.getchild(q::TopDownQuadtreeCursor, i::Int)
    n = STI.nchild(q)
    1 <= i <= n || throw(BoundsError(q, i))
    i_is_leaf = length(q.irange) <= q.leafsize[1]
    j_is_leaf = length(q.jrange) <= q.leafsize[2]
    if i_is_leaf
        s = length(q.jrange) ÷ 2
        return _child(q, q.irange, i == 1 ? _lo(q.jrange, s) : _hi(q.jrange, s))
    elseif j_is_leaf
        s = length(q.irange) ÷ 2
        return _child(q, i == 1 ? _lo(q.irange, s) : _hi(q.irange, s), q.jrange)
    else
        si = length(q.irange) ÷ 2
        sj = length(q.jrange) ÷ 2
        return _child(q, i <= 2 ? _lo(q.irange, si) : _hi(q.irange, si), isodd(i) ? _lo(q.jrange, sj) : _hi(q.jrange, sj))
    end
end
STI.getchild(q::TopDownQuadtreeCursor) = (STI.getchild(q, i) for i in 1:STI.nchild(q))

STI.node_extent(q::TopDownQuadtreeCursor) = blockextent(q.grid, q.face, q.irange, q.jrange)

@inline _emit(q::TopDownQuadtreeCursor{G,false}, i, j) where {G} = faceposition(q.grid, q.face, i, j)
@inline _emit(q::TopDownQuadtreeCursor{G,true}, i, j) where {G} = (i - q.i0 + 1) + (j - q.j0) * q.ni

@inline function _leaf_entry(q::TopDownQuadtreeCursor, k::Int)
    a, b = divrem(k - 1, length(q.irange))
    i = first(q.irange) + b
    j = first(q.jrange) + a
    return (_emit(q, i, j), blockextent(q.grid, q.face, i:i, j:j))
end

"""
    LeafEntries

A leaf's `(position, extent)` pairs, up to `C` of them, built once per visit
with no heap storage (CR's dual DFS carries a leaf's entries while the other
tree descends).
"""
struct LeafEntries{C,E} <: AbstractVector{Tuple{Int,E}}
    n::Int
    e::NTuple{C,Tuple{Int,E}}
end
Base.size(l::LeafEntries) = (l.n,)
Base.IndexStyle(::Type{<:LeafEntries}) = IndexLinear()
Base.@propagate_inbounds function Base.getindex(l::LeafEntries, k::Int)
    @boundscheck 1 <= k <= l.n || throw(BoundsError(l, k))
    return l.e[k]
end

function STI.child_indices_extents(q::TopDownQuadtreeCursor{G,L,C}) where {G,L,C}
    n = length(q.irange) * length(q.jrange)
    e1 = _leaf_entry(q, 1)
    es = ntuple(k -> k == 1 ? e1 : k <= n ? _leaf_entry(q, k) : e1, Val(C))
    return LeafEntries{C,typeof(e1[2])}(n, es)
end

split_weight(q::TopDownQuadtreeCursor) = length(q.irange) * length(q.jrange)

# the leaf walk over the quadrants, emitting positions without extents
function _leaf_walk!(out::Vector{Int}, q::TopDownQuadtreeCursor)
    if STI.isleaf(q)
        for j in q.jrange, i in q.irange
            push!(out, _emit(q, i, j))
        end
    else
        for k in 1:STI.nchild(q)
            _leaf_walk!(out, STI.getchild(q, k))
        end
    end
    return out
end
subtree_indices(q::TopDownQuadtreeCursor) = _leaf_walk!(sizehint!(Int[], split_weight(q)), q)

"""
    block_indices(q::TopDownQuadtreeCursor)
    block_indices(t::FaceRootTree)

The positions under a node as one lazy ascending block (`TileIndices`, or the
tile-local numbering for a tile's own tree; `Base.OneTo(ncells)` at a face
root): the same set as `subtree_indices`, in position order rather than leaf
order.
"""
block_indices(q::TopDownQuadtreeCursor{G,false}) where {G} = TileIndices(q.grid, q.face, q.irange, q.jrange)
block_indices(q::TopDownQuadtreeCursor{G,true}) where {G} =
    LocalTileIndices(q.irange, q.jrange, q.i0, q.j0, q.ni)

"""
    FaceRootTree(grid, leafsize = (2, 2))

The tree of a multi-face structured grid: a root whose children are the
[`TopDownQuadtreeCursor`](@ref)s of the faces. Its extent is the whole sphere
on a curved manifold, the union of the faces' boxes on the plane.
"""
struct FaceRootTree{G}
    grid::G
    leafsize::NTuple{2,Int}
end
FaceRootTree(grid) = FaceRootTree(grid, (2, 2))
STI.isspatialtree(::Type{<:FaceRootTree}) = true
STI.node_extent_is_expensive(::Type{<:FaceRootTree}) = true
STI.isleaf(::FaceRootTree) = false
STI.nchild(t::FaceRootTree) = nfaces(t.grid)
STI.getchild(t::FaceRootTree, k::Int) = TopDownQuadtreeCursor(t.grid, k; leafsize=t.leafsize)
STI.getchild(t::FaceRootTree) = (STI.getchild(t, k) for k in 1:STI.nchild(t))
function STI.node_extent(t::FaceRootTree)
    m = best_manifold(t.grid)
    iscurved(m) && return _whole_sphere_cap(Float64)
    return reduce(Extents.union, (STI.node_extent(c) for c in STI.getchild(t)))
end
split_weight(t::FaceRootTree) = ncells(t.grid)
function subtree_indices(t::FaceRootTree)
    out = sizehint!(Int[], ncells(t.grid))
    for k in 1:STI.nchild(t)
        _leaf_walk!(out, STI.getchild(t, k))
    end
    return out
end
block_indices(t::FaceRootTree) = Base.OneTo(ncells(t.grid))
Base.show(io::IO, t::FaceRootTree) = print(io, "FaceRootTree(", nfaces(t.grid), " faces)")

# A block's span under the face numbering: one run per row when positions
# ascend `i` fastest (a `LatticeBlock`, exact against a range), the local
# numbering for a tile's own tree, unknown for any other data order.
GrI.position_span(q::TopDownQuadtreeCursor{G,false}) where {G} =
    facecolumnmajor(q.grid, q.face) ?
        GrI.LatticeBlock(faceposition(q.grid, q.face, first(q.irange), first(q.jrange)),
            length(q.irange), facesize(q.grid, q.face)[1], length(q.jrange)) :
        nothing
GrI.position_span(q::TopDownQuadtreeCursor{G,true}) where {G} =
    GrI.LatticeBlock(_emit(q, first(q.irange), first(q.jrange)), length(q.irange), q.ni, length(q.jrange))
GrI.position_span(t::FaceRootTree) = 1:ncells(t.grid)

# A restricted leaf keeps the tuple container: members are compacted into the
# same `NTuple`, so a masked quadtree leaf still never reaches the heap.
function GrI._members(es::LeafEntries, cells, ::Val{L}) where {L}
    out = es.e
    n = 0
    for (p, e) in es
        k = GrI._memberindex(cells, p)
        k == 0 && continue
        n += 1
        out = Base.setindex(out, (L ? k : p, e), n)
    end
    return LeafEntries(n, out)
end
