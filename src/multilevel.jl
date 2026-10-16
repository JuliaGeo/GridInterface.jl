# ===========================================================================
# A multi-level cell set: an ordinary grid whose cells come from several levels
# of one system. Parent system, not parent grid.
#
# Two layouts:
#  - sorted (has_sorted_subtrees): cells in curve order at the reference level,
#    keys[i] = first(descendant_range(sys, cells[i], ref)), strictly ascending.
#    localindex is one descendant_range + a binary search; the tree is the
#    window-mode HierarchyCursor; cellat goes through the reference level.
#  - unsorted: cells sorted by (level, id) with a level-offset table; localindex
#    is a binary search in that level's run; the tree is the SelectionCursor.
# ===========================================================================

"""
    MultiLevelSet(sys, cells; check = true)

A grid over cells drawn from several levels of system `sys`. Law: no member is
an ancestor of another (checked in O(n) on sorted systems). Geometry forwards
to `levelgrid(sys, level(c))`; identity is the typed id; it answers `system`
and refuses `globalindex`; its tessellation trait is the system's refinement
trait; it is non-conforming (hanging vertices) so it has no dual.

Two regrid behaviours, both available:

- **direct**: pass the set itself; each member is one column (one polygon).
  New behaviour, 2-7x faster; differs from the expansion where a coarse
  member's polygon differs from the union of its fine descendants
  (non-congruent systems, densified boundaries).
- **expanded**: pass [`expand_to_reference`](@ref)`(set)`, the subset of the
  reference level holding every member's descendants; bit-identical to today's
  expansion path. [`aggregation_pairs`](@ref) maps its columns back to members.
"""
struct MultiLevelSet{S,ID} <: AbstractGrid
    system::S
    cells::Vector{ID}
    keys::Vector{Int}         # sorted layout; empty otherwise
    offsets::Vector{Int}      # unsorted layout: run of level l is offsets[l-l0+1]:offsets[l-l0+2]-1
    reference_level::Int
    sorted::Bool
end

function MultiLevelSet(sys, cells::AbstractVector; check::Bool=true)
    isempty(cells) && throw(ArgumentError("a multi-level set needs at least one cell"))
    ID = eltype(cells)
    ref = maximum(level, cells)
    if has_sorted_subtrees(sys)
        keys = [first(descendant_range(sys, c, ref)) for c in cells]
        perm = sortperm(keys)
        cs = Vector{ID}(cells[perm]); ks = keys[perm]
        s = MultiLevelSet{typeof(sys),ID}(sys, cs, ks, Int[], ref, true)
        check && (check_antichain(s) || throw(ArgumentError("a member is an ancestor of another (or duplicated)")))
        return s
    end
    l0 = first(levels(sys))
    cs = sort(Vector{ID}(cells); by=c -> (level(c), c))
    offsets = zeros(Int, ref - l0 + 2)
    for c in cs
        offsets[level(c)-l0+2] += 1
    end
    offsets[1] = 1
    offsets = cumsum(offsets)
    s = MultiLevelSet{typeof(sys),ID}(sys, cs, Int[], offsets, ref, false)
    check && (check_antichain(s) || throw(ArgumentError("a member is an ancestor of another (or duplicated)")))
    return s
end

"""
    check_antichain(s::MultiLevelSet) -> Bool

No member is an ancestor of another. O(n) from the keys on sorted systems,
O(n · depth) with `ancestor` otherwise.
"""
function check_antichain(s::MultiLevelSet)
    if s.sorted
        for i in 2:length(s.cells)
            prev = descendant_range(s.system, s.cells[i-1], s.reference_level)
            s.keys[i] > last(prev) || return false
        end
        return true
    end
    members = Set(s.cells)
    length(members) == length(s.cells) || return false
    for c in s.cells, l in first(levels(s.system)):level(c)-1
        ancestor(s.system, c, l) in members && return false
    end
    return true
end

system(s::MultiLevelSet) = s.system
gridrelation(::MultiLevelSet) = ParentSystem()
ncells(s::MultiLevelSet) = length(s.cells)
Base.@propagate_inbounds function cellindex(s::MultiLevelSet, i::Integer)
    @boundscheck (1 <= i <= length(s.cells)) || _throw_bounds(s, i)
    return @inbounds s.cells[i]
end
cellindextypes(s::MultiLevelSet{S,ID}) where {S,ID} = (ID,)

@inline _lg(s::MultiLevelSet, c) = levelgrid(s.system, level(c))
_refgrid(s::MultiLevelSet) = levelgrid(s.system, s.reference_level)

GOC.best_manifold(s::MultiLevelSet) = best_manifold(_refgrid(s))
cellshape(s::MultiLevelSet) = cellshape(_refgrid(s))
edgekind(s::MultiLevelSet) = edgekind(_refgrid(s))
cell_geometry_is_expensive(s::MultiLevelSet) = cell_geometry_is_expensive(_refgrid(s))
tessellation(s::MultiLevelSet) = has_congruent_refinement(s.system) ? Tessellating() : Overlapping()
conformity(s::MultiLevelSet) = NonConforming()
# Analytic through the reference level. Under non-congruent refinement members
# overlap and this is the hierarchy's answer (one containing member, not
# necessarily the lowest position); the PRD lets an analytic answer win there.
location_is_analytic(s::MultiLevelSet) = s.sorted && location_is_analytic(_refgrid(s))
extent_is_expensive(s::MultiLevelSet) = extent_is_expensive(_refgrid(s))

# Geometry: one forwarding step, type stable iff `levelgrid(sys, ::Int)` has one
# concrete return type for every level.
for f in (:cell_boundary, :cell_geometry, :cell_centroid, :cell_area, :cell_extent, :cell_corners)
    @eval Base.@propagate_inbounds function $f(s::MultiLevelSet, i::Integer)
        c = cellindex(s, i)
        return $f(_lg(s, c), c)
    end
end
function cell_contains(s::MultiLevelSet, i::Integer, p)
    c = cellindex(s, i)
    g = _lg(s, c)
    return cell_contains(g, localindex(g, c), p)
end

function localindex(s::MultiLevelSet, c)
    level(c) <= s.reference_level || return nothing
    if s.sorted
        k = first(descendant_range(s.system, c, s.reference_level))
        i = searchsortedfirst(s.keys, k)
        (i <= length(s.cells) && @inbounds(s.cells[i]) == c) || return nothing
        return i
    end
    l0 = first(levels(s.system))
    li = level(c) - l0 + 1
    (1 <= li < length(s.offsets)) || return nothing
    lo, hi = s.offsets[li], s.offsets[li+1] - 1
    i = searchsortedfirst(s.cells, c, lo, hi, Base.Order.By(x -> (level(x), x)))
    (i <= hi && @inbounds(s.cells[i]) == c) || return nothing
    return i
end
localindex(s::MultiLevelSet, c::Integer) = throw(ArgumentError("a multi-level set has no position-named ids; pass a typed id"))

"""
    member_at_reference(s, p) -> Union{Int,Nothing}

The member whose subtree holds reference-level position `p`.
"""
@inline function member_at_reference(s::MultiLevelSet, p::Int)
    i = searchsortedlast(s.keys, p)
    i >= 1 || return nothing
    return p <= last(descendant_range(s.system, @inbounds(s.cells[i]), s.reference_level)) ? i : nothing
end

# Analytic location (sorted layout over an analytic reference level): locate
# at the reference level, then map to the member. Exact geometrically only
# under congruent refinement; elsewhere it is the hierarchy's answer, which
# wins over the geometric tie rule (PRD).
function cellat(s::MultiLevelSet, p)
    if location_is_analytic(s)
        q = cellat(_refgrid(s), p)
        q === nothing && return nothing
        return member_at_reference(s, q)
    end
    return fallback_cellat(best_manifold(s), s, p)
end

# --- frontier and tree ------------------------------------------------------
fr_keys(s::MultiLevelSet) = s.keys
fr_system(s::MultiLevelSet) = s.system
fr_ref(s::MultiLevelSet) = s.reference_level
fr_n(s::MultiLevelSet) = length(s.cells)
@inline fr_level(s::MultiLevelSet, k) = level(@inbounds s.cells[k])
@inline fr_id(s::MultiLevelSet, k) = @inbounds s.cells[k]
@inline fr_emit(s::MultiLevelSet, k) = k
@inline fr_extent(s::MultiLevelSet, k) = cell_extent(s, k)

# `has_sorted_subtrees` is a type-level fact for every system seen so far, so
# this branch folds and the return type is concrete.
treeify(s::MultiLevelSet, bucket::Integer) =
    has_sorted_subtrees(s.system) ? HierarchyCursor(s; bucket) : SelectionCursor(s; bucket)

_fr_extent_is_expensive(s::MultiLevelSet) = extent_is_expensive(s)

"""
    MLSFrontier(set, positions, emit)

A window of a sorted multi-level set (ascending set positions): the frontier
behind `treeify(subset(set, …))` (emits local positions) and
`subtree(set, cells)` (emits set positions). Window mode, O(1) to build.
"""
struct MLSFrontier{M,V<:AbstractVector{Int},B}
    s::M
    pos::V
end
MLSFrontier(s, pos, emit::Bool) = MLSFrontier{typeof(s),typeof(pos),emit}(s, pos)
fr_keys(F::MLSFrontier) = view(F.s.keys, F.pos)
fr_system(F::MLSFrontier) = F.s.system
fr_ref(F::MLSFrontier) = F.s.reference_level
fr_n(F::MLSFrontier) = length(F.pos)
@inline fr_level(F::MLSFrontier, k) = level(@inbounds F.s.cells[F.pos[k]])
@inline fr_id(F::MLSFrontier, k) = @inbounds F.s.cells[F.pos[k]]
@inline fr_emit(F::MLSFrontier{M,V,true}, k) where {M,V} = @inbounds F.pos[k]
@inline fr_emit(F::MLSFrontier{M,V,false}, k) where {M,V} = k
@inline fr_emitted(F::MLSFrontier{M,V,true}, ks::AbstractVector{Int}) where {M,V} = view(F.pos, ks)
@inline fr_extent(F::MLSFrontier, k) = cell_extent(F.s, @inbounds F.pos[k])
_fr_extent_is_expensive(F::MLSFrontier) = extent_is_expensive(F.s)

function _mls_cursor(m::MultiLevelSet, pos, emit::Bool; root=nothing, bucket=1)
    m.sorted || return fallback_subtree(best_manifold(m), m, pos, bucket)
    return HierarchyCursor(MLSFrontier(m, pos, emit); root, bucket)
end

treeify(s::Subset{<:MultiLevelSet}, n::Integer) =
    s.parent.sorted ? _mls_cursor(s.parent, s.cells, false; root=s.root, bucket=n) : _default_treeify(s, Int(n))
subtree(m::MultiLevelSet, cells::AbstractVector{<:Integer}) =
    _mls_cursor(m, cells, true; bucket=leaf_capacity(m, Balanced()))
subtree(m::MultiLevelSet, s::Subset) = parent(s) === m ?
    _mls_cursor(m, s.cells, true; root=s.root, bucket=leaf_capacity(m, Balanced())) :
    subtree(m, parentpositions(s))

"""
    chunks(set::MultiLevelSet, l) -> AbstractVector of subsets

The members grouped by their level-`l` ancestor (one window of the sorted
layout per ancestor, rooted at it); a member at or above level `l` is its own
chunk, rooted at itself. Sorted systems only.
"""
function hierarchy_chunks(m::MultiLevelSet, l::Integer)
    m.sorted || throw(ArgumentError("chunks of a multi-level set need sorted subtrees"))
    sys, ref = m.system, m.reference_level
    ID = eltype(m.cells)
    starts = Int[]; roots = ID[]
    n = length(m.cells)
    k = 1
    while k <= n
        c = m.cells[k]
        push!(starts, k)
        if level(c) <= l
            push!(roots, c)
            k += 1
        else
            a = ancestor(sys, c, l)
            r = descendant_range(sys, a, ref)
            push!(roots, a)
            k = searchsortedlast(m.keys, last(r), k, n, Base.Order.Forward) + 1
        end
    end
    push!(starts, n + 1)
    return WindowChunks{typeof(m),ID}(m, starts, roots)
end

"""
    expand_to_reference(set::MultiLevelSet) -> Subset

The **expanded** regrid behaviour: the subset of `levelgrid(system, reference
level)` holding every member's descendants, in ascending position. Its
positions are a lazy concatenation of `descendant_range`s (storage O(members));
its tree is the hierarchy cursor. Sorted systems only.
"""
function expand_to_reference(m::MultiLevelSet)
    m.sorted || throw(ArgumentError("expand_to_reference needs sorted subtrees"))
    rs = UnitRange{Int}[descendant_range(m.system, c, m.reference_level) for c in m.cells]
    return Subset(_refgrid(m), RangeCat(rs))
end

"""
    aggregation_pairs(s::MultiLevelSet) -> (rows, cols)

The 0/1 expansion of members to their reference-level descendants, as
`(rows, cols)` index vectors: row = reference-level position, col = member.
(Size of the expansion, by definition.)
"""
function aggregation_pairs(s::MultiLevelSet)
    s.sorted || throw(ArgumentError("needs sorted subtrees"))
    rows = Int[]; cols = Int[]
    for (i, c) in pairs(s.cells), r in descendant_range(s.system, c, s.reference_level)
        push!(rows, r); push!(cols, i)
    end
    return rows, cols
end
