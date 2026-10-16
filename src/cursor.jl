# ===========================================================================
# The generic hierarchy cursor: one tree for complete levels, sorted subsets of
# a level, and multi-level sets, written only against the hierarchy primitives.
# It is the *default* `treeify` for level grids; a provider's own tree wins.
#
# A "frontier" supplies, for window index k in 1:n:
#   fr_keys(F)[k]  ascending position at the reference level of the first
#                  descendant of member k (window mode only)
#   fr_level(F,k)  the member's level
#   fr_id(F,k)     the member's typed id
#   fr_emit(F,k)   the position the leaf emits for member k
#   fr_emitted(F,ks) the positions emitted for window indices `ks`, lazily
#   fr_extent(F,k) the member's own extent
# ===========================================================================

# Identity-emitting frontiers answer the window itself; key-emitting ones a
# view of their keys. `subtree_indices` of a cursor is O(1) through this.
@inline fr_emitted(F, ks::AbstractVector{Int}) = ks

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
_fr_extent_is_expensive(F::LevelFrontier) = extent_is_expensive(F.grid)

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
@inline fr_emitted(F::SubsetFrontier{G,S,K,true}, ks::AbstractVector{Int}) where {G,S,K} = view(F.keys, ks)
@inline fr_id(F::SubsetFrontier, k) = cellindex(F.grid, fr_emit(F, k))
@inline fr_extent(F::SubsetFrontier, k) = cell_extent(F.grid, fr_emit(F, k))
_fr_extent_is_expensive(F::SubsetFrontier) = extent_is_expensive(F.grid)

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

# --- the root index: a coarse frontier above a system's roots -----------------

"""
    RootIndex(sys; fanout = 16)

A packed `fanout`-ary index over a system's root cells, in a spatial (STR)
order, with precomputed union extents per group. The generic cursor descends it
between the synthetic root and the roots, so a system with tens of thousands of
roots (CopernicusDEM: 64,800) is not a flat fan-out. Built once, O(R) node
extents; a system returns it from [`root_index`](@ref) (caching it).
"""
struct RootIndex{ID,E}
    ids::Vector{ID}               # root ids in index order
    fanout::Int
    exts::Vector{Vector{E}}       # exts[L][k]: extent of group k at index level L (L = 1 groups roots)
    top::E                        # extent of everything
end

function RootIndex(sys; fanout::Integer=16)
    f = Int(fanout)
    f >= 2 || throw(ArgumentError("fanout must be at least 2"))
    rc = rootcells(sys)
    E = typeof(node_extent(sys, first(rc)))
    ext0 = E[node_extent(sys, c) for c in rc]
    order = _str_order(ext0, f)
    ids = collect(rc)[order]
    cur = ext0[order]
    exts = Vector{E}[]
    while length(cur) > f
        nxt = E[reduce(Extents.union, @view cur[j:min(j + f - 1, end)]) for j in 1:f:length(cur)]
        push!(exts, nxt)
        cur = nxt
    end
    return RootIndex{eltype(ids),E}(ids, f, exts, reduce(Extents.union, cur))
end
_ngroups(ri::RootIndex, L) = L == 0 ? length(ri.ids) : length(ri.exts[L])

# sort-tile-recursive order on extent centres (longitude/z on the sphere)
_centre2(e::SphericalCap) = (atan(Float64(e.point[2]), Float64(e.point[1])), Float64(e.point[3]))
_centre2(e::Extents.Extent) = ((e.X[1] + e.X[2]) / 2, (e.Y[1] + e.Y[2]) / 2)
function _str_order(exts, f)
    n = length(exts)
    c = map(_centre2, exts)
    p = sortperm(c; by=first)
    nslab = max(1, ceil(Int, sqrt(n / f)))
    slab = cld(n, nslab)
    for j in 1:slab:n
        r = j:min(j + slab - 1, n)
        p[r] = p[r][sortperm(c[p[r]]; by=last)]
    end
    return p
end

"""
    root_index(sys) -> Union{Nothing, RootIndex}

Root-index hook of the generic cursor. `nothing` (default): the synthetic root
fans out to `rootcells(sys)` directly (fine for 12 or 122 roots). A system with
many roots returns a cached [`RootIndex`](@ref) (or a provider-built one).
"""
root_index(sys) = nothing

# --- the window-mode cursor ------------------------------------------------

"""
    HierarchyCursor(frontier; root = nothing, bucket = 1)

STI tree descending a system's own hierarchy from `rootcells` (through the
system's [`root_index`](@ref) when it has one), or from `root`. A node is a
system cell with the contiguous window `lo:hi` of members beneath it, found by
two binary searches over `descendant_range`. It stops at a cell when its window
holds `<= bucket` members or exactly that cell. Requires `has_sorted_subtrees`.
Building it is O(1).
"""
struct HierarchyCursor{F,ID,E,RI}
    fr::F
    top::Int        # first(levels(sys)); a node at level < top is synthetic (root or group)
    level::Int
    id::ID
    lo::Int
    hi::Int
    bucket::Int
    ri::RI          # nothing, or the system's RootIndex
    gl::Int         # group level in `ri` (0: a cell node or a flat synthetic root)
    gk::Int         # group number at that level
end

function HierarchyCursor(fr::F; root=nothing, bucket::Integer=1) where {F}
    sys = fr_system(fr)
    top = first(levels(sys))
    E = Base.promote_op(fr_extent, F, Int)
    ri = root === nothing ? root_index(sys) : nothing
    if root === nothing
        id = first(rootcells(sys))
        gl = ri === nothing ? 0 : length(ri.exts) + 1
        return HierarchyCursor{F,typeof(id),E,typeof(ri)}(fr, top, top - 1, id, 1, fr_n(fr), Int(bucket), ri, gl, 1)
    end
    return HierarchyCursor{F,typeof(root),E,typeof(ri)}(fr, top, level(root), root, 1, fr_n(fr), Int(bucket), ri, 0, 0)
end

@inline _count(n::HierarchyCursor) = n.hi - n.lo + 1
@inline _synthetic(n) = n.level < n.top
@inline _isgroup(n::HierarchyCursor) = n.ri !== nothing && n.gl >= 1

STI.isspatialtree(::Type{<:HierarchyCursor}) = true
STI.node_extent_is_expensive(::Type{<:HierarchyCursor}) = true
function split_weight(n::HierarchyCursor)
    c = max(0, _count(n))
    _isgroup(n) || return c
    ri = n.ri
    nr = length(ri.ids)
    span = min(ri.fanout^min(n.gl, 8), nr)     # roots under one group at this level
    return max(1, (c * span) ÷ nr)
end
# A root-index group keeps its parent's window, so its leaves are gathered from
# the windows of its own roots, in index order (roots are in STR order).
function subtree_indices(n::HierarchyCursor)
    _isgroup(n) || return fr_emitted(n.fr, n.lo:n.hi)
    ri = n.ri
    nr = length(ri.ids)
    span = n.gl > length(ri.exts) ? nr : ri.fanout^n.gl   # roots under one group at this level
    out = Int[]
    for k in (n.gk - 1) * span + 1:min(n.gk * span, nr)
        c = _child(n, @inbounds ri.ids[k])
        append!(out, fr_emitted(n.fr, c.lo:c.hi))
    end
    return out
end

function STI.isleaf(n::HierarchyCursor)
    k = _count(n)
    k <= 0 && return true
    _isgroup(n) && return n.gl <= length(n.ri.exts) ? false : k <= n.bucket
    k <= n.bucket && return true
    n.level >= fr_ref(n.fr) && return true
    return k == 1 && fr_level(n.fr, n.lo) <= n.level
end

# The ids children are taken from: rootcells at a flat synthetic root, the root
# index's ids below a group, the system's children otherwise.
function _child_ids(n::HierarchyCursor)
    _isgroup(n) && return n.ri.ids
    return _synthetic(n) ? rootcells(fr_system(n.fr)) : children(fr_system(n.fr), n.id)
end
function _child_range(n::HierarchyCursor, ids)
    _isgroup(n) || return 1:length(ids)
    f = n.ri.fanout
    nk = _ngroups(n.ri, n.gl - 1)
    lo = n.gl > length(n.ri.exts) ? 1 : (n.gk - 1) * f + 1
    hi = n.gl > length(n.ri.exts) ? nk : min(n.gk * f, nk)
    return lo:hi
end
@inline function _kth_child(n::HierarchyCursor{F,ID,E,RI}, ids, k) where {F,ID,E,RI}
    if _isgroup(n) && n.gl >= 2
        return HierarchyCursor{F,ID,E,RI}(n.fr, n.top, n.level, n.id, n.lo, n.hi, n.bucket, n.ri, n.gl - 1, k)
    end
    return _child(n, @inbounds ids[k])
end

@inline function _child(n::HierarchyCursor{F,ID,E,RI}, cid) where {F,ID,E,RI}
    fr = n.fr
    r = descendant_range(fr_system(fr), cid, fr_ref(fr))
    ks = fr_keys(fr)
    lo = searchsortedfirst(ks, first(r), n.lo, n.hi, Base.Order.Forward)
    hi = searchsortedlast(ks, last(r), n.lo, n.hi, Base.Order.Forward)
    return HierarchyCursor{F,ID,E,RI}(fr, n.top, n.level + 1, cid, lo, hi, n.bucket, n.ri, 0, 0)
end
_nonempty(n) = _count(n) > 0

function STI.getchild(n::HierarchyCursor)
    ids = _child_ids(n)
    return Iterators.filter(_nonempty, (_kth_child(n, ids, k) for k in _child_range(n, ids)))
end
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
    if _isgroup(n)
        ri = n.ri
        return n.gl > length(ri.exts) ? ri.top : @inbounds ri.exts[n.gl][n.gk]
    end
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

# Cheap extents stay lazy (no allocation). Expensive ones are built once per
# visit: CR's dual DFS carries a leaf's entries while the opposing tree
# descends, and a lazy vector would recompute every cap per opposing leaf.
function STI.child_indices_extents(n::HierarchyCursor{F,ID,E}) where {F,ID,E}
    entries = LeafEntries{F,E}(n.fr, n.lo, max(0, _count(n)))
    return _fr_extent_is_expensive(n.fr) ? collect(entries) : entries
end

# --- the selection-mode cursor (systems without sorted subtrees) -------------

"""
    SelectionCursor(frontier; root = nothing, bucket = 1)

The variant for systems without sorted subtrees: a node carries the window
indices of its members, and children are found by bucketing them with
`ancestor`. O(n · depth) ancestor calls and one vector per node: the one tree
in the package that is not O(1) to build.
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
subtree_indices(n::SelectionCursor) = fr_emitted(n.fr, n.sel)

function STI.isleaf(n::SelectionCursor)
    k = length(n.sel)
    k <= n.bucket && return true
    n.level >= fr_ref(n.fr) && return true
    return k == 1 && fr_level(n.fr, n.sel[1]) <= n.level
end

_child_ids(n::SelectionCursor) = _synthetic(n) ? rootcells(fr_system(n.fr)) : children(fr_system(n.fr), n.id)

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
    hierarchy_subtree(grid, cells; bucket = 1, root = nothing)

The hierarchy cursor over positions `cells` (ascending) of a complete level
grid, emitting those parent positions; `root` starts it at a known ancestor.
"""
function hierarchy_subtree(grid, cells::AbstractVector{<:Integer}; bucket::Integer=1, root=nothing)
    sys = system(grid)
    keys = cells isa AbstractVector{Int} ? cells : Int.(cells)
    (keys isa AbstractRange || keys isa RangeCat || issorted(keys)) || (keys = sort(keys))
    fr = SubsetFrontier(grid, sys, level(grid), keys, true)
    return _cursor(fr, sys; root, bucket)
end

# `subtree` dispatches on the collection: a chunk (a Subset of this very grid)
# with a root starts at that root, O(1).
function subtree(grid, s::Subset)
    if parent(s) === grid && s.root !== nothing && system(grid) !== nothing && level(grid) !== nothing &&
       gridrelation(grid) isa NoParent
        return hierarchy_subtree(grid, s.cells; root=s.root, bucket=_bucket(grid))
    end
    return subtree(grid, parentpositions(s))
end

treeify(s::Subset) = _default_treeify(s)
