# ===========================================================================
# The generic hierarchy cursor: one tree for complete levels, sorted subsets of
# a level, and multi-level sets, written only against the hierarchy primitives.
#
# A "frontier" supplies, for window index k in 1:n:
#   fr_key(F)[k]   ascending position at the reference level of the first
#                  descendant of member k (window mode only)
#   fr_level(F,k)  the member's level
#   fr_id(F,k)     the member's typed id (selection mode only)
#   fr_emit(F,k)   the position the leaf emits for member k
#   fr_extent(F,k) the member's own extent
# ===========================================================================

"""
    LevelFrontier(grid, sys, ref)

Members are every cell of a complete level grid; key == position.
"""
struct LevelFrontier{G,S}
    grid::G
    sys::S
    ref::Int
    n::Int
end
LevelFrontier(grid, sys, ref::Integer) = LevelFrontier(grid, sys, Int(ref), ncells(grid))
fr_keys(F::LevelFrontier) = Base.OneTo(F.n)
@inline fr_level(F::LevelFrontier, k) = F.ref
@inline fr_id(F::LevelFrontier, k) = cellindex(F.grid, k)
@inline fr_emit(F::LevelFrontier, k) = k
@inline fr_extent(F::LevelFrontier, k) = cell_extent(F.grid, k)

"""
    SubsetFrontier(grid, sys, ref, keys, emitkeys)

Members are cells of one level at reference-level positions `keys` (ascending).
Leaves emit `keys[k]` when `emitkeys` (a `subtree` of a level grid: parent
positions), else `k` (a subset grid: its own positions). Extents are asked of
`grid` at the emitted position.
"""
struct SubsetFrontier{G,S,K<:AbstractVector{Int},B}
    grid::G
    sys::S
    ref::Int
    keys::K
end
SubsetFrontier(grid, sys, ref, keys, emitkeys::Bool) =
    SubsetFrontier{typeof(grid),typeof(sys),typeof(keys),emitkeys}(grid, sys, Int(ref), keys)
fr_keys(F::SubsetFrontier) = F.keys
@inline fr_level(F::SubsetFrontier, k) = F.ref
@inline fr_emit(F::SubsetFrontier{G,S,K,true}, k) where {G,S,K} = @inbounds F.keys[k]
@inline fr_emit(F::SubsetFrontier{G,S,K,false}, k) where {G,S,K} = k
@inline fr_id(F::SubsetFrontier, k) = cellindex(F.grid, fr_emit(F, k))
@inline fr_extent(F::SubsetFrontier, k) = cell_extent(F.grid, fr_emit(F, k))

fr_n(F) = length(fr_keys(F))
fr_system(F) = F.sys
fr_ref(F) = F.ref

"""
    GlobalKeys(grid)

Lazy `globalindex(grid, cellindex(grid, k))` over `1:ncells(grid)`: the keys of
any subset type that answers `globalindex`, without materialising them.
"""
struct GlobalKeys{G} <: AbstractVector{Int}
    grid::G
    n::Int
end
GlobalKeys(g) = GlobalKeys(g, ncells(g))
Base.size(k::GlobalKeys) = (k.n,)
Base.IndexStyle(::Type{<:GlobalKeys}) = IndexLinear()
Base.getindex(k::GlobalKeys, i::Int) = globalindex(k.grid, cellindex(k.grid, i))::Int

# --- the window-mode cursor ------------------------------------------------

"""
    HierarchyCursor(frontier; root = nothing, bucket = 1)

STI tree descending a system's own hierarchy from `rootcells` (or from `root`).
A node is a system cell with the contiguous window `lo:hi` of members beneath
it, found by two binary searches over `descendant_range`. It stops at a cell
when its window holds `<= bucket` members or exactly that cell. Requires
`has_sorted_subtrees`. Building it is O(1).
"""
struct HierarchyCursor{F,ID,E}
    fr::F
    top::Int        # first(levels(sys)); a node at level < top is the synthetic root
    level::Int
    id::ID
    lo::Int
    hi::Int
    bucket::Int
end

function HierarchyCursor(fr::F; root=nothing, bucket::Integer=1) where {F}
    sys = fr_system(fr)
    top = first(levels(sys))
    E = Base.promote_op(fr_extent, F, Int)
    if root === nothing
        id = first(rootcells(sys))
        return HierarchyCursor{F,typeof(id),E}(fr, top, top - 1, id, 1, fr_n(fr), Int(bucket))
    end
    return HierarchyCursor{F,typeof(root),E}(fr, top, level(root), root, 1, fr_n(fr), Int(bucket))
end

@inline _count(n::HierarchyCursor) = n.hi - n.lo + 1
@inline _synthetic(n) = n.level < n.top

STI.isspatialtree(::Type{<:HierarchyCursor}) = true
STI.node_extent_is_expensive(::Type{<:HierarchyCursor}) = true
split_weight(n::HierarchyCursor) = max(0, _count(n))

function STI.isleaf(n::HierarchyCursor)
    k = _count(n)
    k <= 0 && return true
    k <= n.bucket && return true
    n.level >= fr_ref(n.fr) && return true
    return k == 1 && fr_level(n.fr, n.lo) <= n.level
end

_child_ids(n) = _synthetic(n) ? rootcells(fr_system(n.fr)) : children(fr_system(n.fr), n.id)

@inline function _child(n::HierarchyCursor{F,ID,E}, cid) where {F,ID,E}
    fr = n.fr
    r = descendant_range(fr_system(fr), cid, fr_ref(fr))
    ks = fr_keys(fr)
    lo = searchsortedfirst(ks, first(r), n.lo, n.hi, Base.Order.Forward)
    hi = searchsortedlast(ks, last(r), n.lo, n.hi, Base.Order.Forward)
    return HierarchyCursor{F,ID,E}(fr, n.top, n.level + 1, cid, lo, hi, n.bucket)
end
_nonempty(n) = _count(n) > 0

STI.getchild(n::HierarchyCursor) = Iterators.filter(_nonempty, (_child(n, c) for c in _child_ids(n)))
STI.nchild(n::HierarchyCursor) = count(_ -> true, STI.getchild(n))
function STI.getchild(n::HierarchyCursor, i::Int)
    k = 0
    for c in STI.getchild(n)
        k += 1
        k == i && return c
    end
    throw(BoundsError(n, i))
end

_root_extent(sys) = mapreduce(r -> node_extent(sys, r), Extents.union, rootcells(sys))

function STI.node_extent(n::HierarchyCursor)
    sys = fr_system(n.fr)
    _synthetic(n) && return _root_extent(sys)
    if _count(n) == 1 && fr_level(n.fr, n.lo) == n.level
        return fr_extent(n.fr, n.lo)          # the node is the member: its own extent
    end
    return node_extent(sys, n.id)             # covering law
end

"""
    LeafEntries

Lazy `(emitted position, extent)` pairs over a leaf's window: an
`AbstractVector` with no storage.
"""
struct LeafEntries{F,E} <: AbstractVector{Tuple{Int,E}}
    fr::F
    lo::Int
    n::Int
end
Base.size(e::LeafEntries) = (e.n,)
Base.IndexStyle(::Type{<:LeafEntries}) = IndexLinear()
Base.@propagate_inbounds function Base.getindex(e::LeafEntries, i::Int)
    k = e.lo + i - 1
    return (fr_emit(e.fr, k), fr_extent(e.fr, k))
end

function STI.child_indices_extents(n::HierarchyCursor{F,ID,E}) where {F,ID,E}
    entries = LeafEntries{F,E}(n.fr, n.lo, max(0, _count(n)))
    # cheap extents stay lazy (no allocation); expensive ones are built once per visit
    return _fr_extent_is_expensive(n.fr) ? collect(entries) : entries
end
_fr_extent_is_expensive(F) = extent_is_expensive(F.grid)

# --- the selection-mode cursor (systems without sorted subtrees) -------------

"""
    SelectionCursor(frontier; root = nothing, bucket = 1)

The variant for systems without sorted subtrees: a node carries the window
indices of its members, and children are found by bucketing them with
`ancestor`. O(n · depth) ancestor calls and one vector per node.
"""
struct SelectionCursor{F,ID}
    fr::F
    top::Int
    level::Int
    id::ID
    sel::Vector{Int}
    bucket::Int
end

function SelectionCursor(fr::F; root=nothing, bucket::Integer=1) where {F}
    sys = fr_system(fr)
    top = first(levels(sys))
    sel = collect(1:fr_n(fr))
    if root === nothing
        id = first(rootcells(sys))
        return SelectionCursor{F,typeof(id)}(fr, top, top - 1, id, sel, Int(bucket))
    end
    return SelectionCursor{F,typeof(root)}(fr, top, level(root), root, sel, Int(bucket))
end

STI.isspatialtree(::Type{<:SelectionCursor}) = true
STI.node_extent_is_expensive(::Type{<:SelectionCursor}) = true
split_weight(n::SelectionCursor) = length(n.sel)

function STI.isleaf(n::SelectionCursor)
    k = length(n.sel)
    k <= n.bucket && return true
    n.level >= fr_ref(n.fr) && return true
    return k == 1 && fr_level(n.fr, n.sel[1]) <= n.level
end

function _sel_children(n::SelectionCursor{F,ID}) where {F,ID}
    fr = n.fr
    sys = fr_system(fr)
    ids = _child_ids(n)
    buckets = [Int[] for _ in ids]
    cl = n.level + 1
    for k in n.sel
        lk = fr_level(fr, k)
        lk < cl && continue
        a = lk == cl ? fr_id(fr, k) : ancestor(sys, fr_id(fr, k), cl)
        j = findfirst(==(a), ids)
        j === nothing && throw(ArgumentError("member $(fr_id(fr, k)) is not under node $(n.id)"))
        push!(buckets[j], k)
    end
    return SelectionCursor{F,ID}[SelectionCursor{F,ID}(fr, n.top, cl, ids[j], buckets[j], n.bucket)
                                 for j in eachindex(ids) if !isempty(buckets[j])]
end
STI.getchild(n::SelectionCursor) = _sel_children(n)
STI.getchild(n::SelectionCursor, i::Int) = _sel_children(n)[i]
STI.nchild(n::SelectionCursor) = length(_sel_children(n))

function STI.node_extent(n::SelectionCursor)
    sys = fr_system(n.fr)
    _synthetic(n) && return _root_extent(sys)
    if length(n.sel) == 1 && fr_level(n.fr, n.sel[1]) == n.level
        return fr_extent(n.fr, n.sel[1])
    end
    return node_extent(sys, n.id)
end

STI.child_indices_extents(n::SelectionCursor) = [(fr_emit(n.fr, k), fr_extent(n.fr, k)) for k in n.sel]

# --- entry points used by the generic `treeify` / `subtree` ------------------

_cursor(fr, sys; root=nothing, bucket=1) =
    has_sorted_subtrees(sys) ? HierarchyCursor(fr; root, bucket) : SelectionCursor(fr; root, bucket)

_subset_keys(g) = GlobalKeys(g)
_subset_keys(s::Subset) = gridrelation(s.parent) isa NoParent ? s.cells : GlobalKeys(s)
_subset_root(g) = nothing
_subset_root(s::Subset) = s.root

"""
    hierarchy_tree(grid; bucket = 1)

The generic hierarchy cursor over a complete level grid or a subset of one.
"""
function hierarchy_tree(grid; bucket::Integer=1)
    sys = system(grid)
    ref = level(grid)
    if gridrelation(grid) isa NoParent
        return _cursor(LevelFrontier(grid, sys, ref), sys; bucket)
    end
    fr = SubsetFrontier(grid, sys, ref, _subset_keys(grid), false)
    return _cursor(fr, sys; root=_subset_root(grid), bucket)
end

"""
    hierarchy_subtree(grid, cells; bucket = 1)

The hierarchy cursor over positions `cells` (ascending) of a complete level
grid, emitting those parent positions.
"""
function hierarchy_subtree(grid, cells::AbstractVector{<:Integer}; bucket::Integer=1, root=nothing)
    sys = system(grid)
    keys = cells isa AbstractVector{Int} ? cells : Int.(cells)
    issorted(keys) || (keys = sort(keys))
    fr = SubsetFrontier(grid, sys, level(grid), keys, true)
    return _cursor(fr, sys; root, bucket)
end

# `subtree` of a chunk: start at the chunk's root cell rather than the synthetic root.
subtree(grid, s::Subset) = parent(s) === grid && s.root !== nothing && gridrelation(grid) isa NoParent ?
    hierarchy_subtree(grid, s.cells; root=s.root) : subtree(grid, parentpositions(s))

treeify(s::Subset) = _default_treeify(s)
