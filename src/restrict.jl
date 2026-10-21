# ===========================================================================
# Restriction: a grid's own tree over a subset of its cells. What the generic
# `subtree(grid, cells)` and `treeify(::Subset)` return: the wrapped tree's
# nodes and extents (they still cover), children whose positions cannot reach a
# member dropped, leaves emitting members only.
# ===========================================================================

"""
    position_span(node) -> UnitRange{Int}, a span bound, or nothing

A bound on the parent positions beneath `node`: their lowest-to-highest range,
`nothing` when the tree cannot say cheaply (an R-tree), or a tighter *span
bound* for a node whose positions are not one run: a [`LatticeBlock`](@ref)
for a rectangle of a column-major lattice, or any value answering
[`span_range`](@ref), [`span_mayhold`](@ref) and [`span_inside`](@ref). A
bound, never the exact set: [`MaskedCursor`](@ref) drops a child only when
its bound misses every member, so a loose bound costs visits, never
correctness. Default `nothing`; every tree GridInterface builds answers.
"""
position_span(node) = nothing

"""
    span_range(bound) -> UnitRange{Int}

The lowest-to-highest range of a span bound (identity on a range).
"""
span_range(r::UnitRange{Int}) = r

"""
    span_mayhold(bound, cells) -> Bool

Whether a member of the strictly ascending `cells` may lie under a span
bound; `true` for `nothing`. Default: the range's overlap with `cells`, O(1)
against a range, one binary search against a vector. A bound answers tighter
where it can (a [`LatticeBlock`](@ref) against a range).
"""
span_mayhold(::Nothing, cells) = true
span_mayhold(b, cells) = _mayhold(span_range(b), cells)
_mayhold(sp::UnitRange{Int}, cells::AbstractUnitRange) =
    max(first(sp), first(cells)) <= min(last(sp), last(cells))
function _mayhold(sp::UnitRange{Int}, cells::AbstractVector)
    isempty(sp) && return false
    k = searchsortedfirst(cells, first(sp))
    return k <= length(cells) && @inbounds(cells[k]) <= last(sp)
end

"""
    span_inside(bound, cells) -> Bool

Whether every position under a span bound is a member of `cells`: only a
range of members can say in O(1); `false` for `nothing` or a vector.
"""
span_inside(::Nothing, cells) = false
span_inside(b, cells) = _inside(span_range(b), cells)
_inside(sp::UnitRange{Int}, cells::AbstractUnitRange) = first(cells) <= first(sp) && last(sp) <= last(cells)
_inside(sp, cells) = false

"""
    LatticeBlock(first, width, stride, rows)

The span bound of a rectangle of a column-major lattice: `rows` runs of
`width` positions, run `k` (from 0) starting at `first + k * stride`. Exact
against a range of members: a block's lowest-to-highest range would cover the
other columns of its rows, and a tree that bisects columns before rows (every
quadtree over a wide face) would then never prune.
"""
struct LatticeBlock
    first::Int
    width::Int
    stride::Int
    rows::Int
end
span_range(b::LatticeBlock) = b.first:(b.first + (b.rows - 1) * b.stride + b.width - 1)
# the runs that end at or after the first member and start at or before the last
function span_mayhold(b::LatticeBlock, cells::AbstractUnitRange)
    (b.rows > 0 && b.width > 0 && !isempty(cells)) || return false
    klo = max(0, cld(first(cells) - (b.first + b.width - 1), b.stride))
    khi = min(b.rows - 1, fld(last(cells) - b.first, b.stride))
    return klo <= khi
end
# Against a vector, exact too (its lowest-to-highest range would keep every
# column half of a wide block, so a sparse restriction would never prune nor
# root below the lattice): from the first member on, either it lies in a run
# or the search jumps to the next run's start. O(min(rows, members) log n),
# and it stops at the first hit.
function span_mayhold(b::LatticeBlock, cells::AbstractVector)
    (b.rows > 0 && b.width > 0) || return false
    n = length(cells)
    k = searchsortedfirst(cells, b.first)
    while k <= n
        p = @inbounds cells[k]
        r = fld(p - b.first, b.stride)            # the run at or before `p`
        r < b.rows || return false
        p - (b.first + r * b.stride) < b.width && return true
        r + 1 < b.rows || return false
        k = searchsortedfirst(cells, b.first + (r + 1) * b.stride, k + 1, n, Base.Order.Forward)
    end
    return false
end

"""
    MaskedCursor{L}(tree, cells)

`tree` (any SpatialTreeInterface node) restricted to the strictly ascending
parent positions `cells`: the same nodes and extents, children whose
[`position_span`](@ref) misses every member dropped, leaves emitting only
members, as the parent positions when `L` is `false` (a subtree in the grid)
or as their index in `cells` when `true` (a subset's own positions). The
constructor descends from `tree` while exactly one child may hold members, so
a chunk roots at its own node in O(depth). `split_weight` is the span's
overlap with `cells` when both are ranges, else the wrapped node's;
`subtree_indices` the members under the node.
"""
struct MaskedCursor{Q,V<:AbstractVector{Int},L}
    node::Q
    cells::V
end
function MaskedCursor{L}(tree, cells::V) where {V<:AbstractVector{Int},L}
    node = _tight_root(tree, cells)
    return MaskedCursor{typeof(node),V,L}(node, cells)
end
_masked(m::MaskedCursor{Q,V,L}, node) where {Q,V,L} = MaskedCursor{typeof(node),V,L}(node, m.cells)

Base.show(io::IO, m::MaskedCursor) = print(io, "MaskedCursor(", m.node, ", ", length(m.cells), " members)")

function _tight_root(node, cells)
    while !STI.isleaf(node)
        only = nothing
        for c in STI.getchild(node)
            span_mayhold(position_span(c), cells) || continue
            only === nothing || return node       # two candidates: this is the root
            only = c
        end
        only === nothing && return node           # no member below: stay, emit nothing
        node = only
    end
    return node
end

# a member's index in `cells`, or 0
@inline function _memberindex(cells, p::Int)
    k = searchsortedfirst(cells, p)
    return (k <= length(cells) && @inbounds(cells[k]) == p) ? k : 0
end
@inline _memberindex(cells::AbstractUnitRange, p::Int) = p in cells ? p - first(cells) + 1 : 0

STI.isspatialtree(::Type{<:MaskedCursor}) = true
STI.node_extent_is_expensive(::Type{<:MaskedCursor{Q}}) where {Q} = STI.node_extent_is_expensive(Q)
STI.isleaf(m::MaskedCursor) = STI.isleaf(m.node)
STI.node_extent(m::MaskedCursor) = STI.node_extent(m.node)
_keeps(m::MaskedCursor, c) = span_mayhold(position_span(c), m.cells)
STI.getchild(m::MaskedCursor) = (_masked(m, c) for c in STI.getchild(m.node) if _keeps(m, c))
STI.nchild(m::MaskedCursor) = count(c -> _keeps(m, c), STI.getchild(m.node))
function STI.getchild(m::MaskedCursor, k::Int)
    i = 0
    for c in STI.getchild(m)
        i += 1
        i == k && return c
    end
    throw(BoundsError(m, k))
end

position_span(m::MaskedCursor{Q,V,false}) where {Q,V} = position_span(m.node)
function position_span(m::MaskedCursor{Q,V,true}) where {Q,V}
    sp = position_span(m.node)
    sp === nothing && return nothing
    return _localspan(span_range(sp), m.cells)
end
# the span mapped to indices in `cells`: still a bound, exact at the ends
_localspan(sp, cells::AbstractUnitRange) =
    (max(first(sp), first(cells)) - first(cells) + 1):(min(last(sp), last(cells)) - first(cells) + 1)
_localspan(sp, cells) = searchsortedfirst(cells, first(sp)):searchsortedlast(cells, last(sp))

function split_weight(m::MaskedCursor)
    sp = position_span(m.node)
    (sp !== nothing && m.cells isa AbstractUnitRange) && return length(intersect(span_range(sp), m.cells))
    return split_weight(m.node)
end

# A leaf whose span lies inside a range of members needs no filter: a subtree
# hands the wrapped entries back as they are, a subset shifts them.
function STI.child_indices_extents(m::MaskedCursor{Q,V,L}) where {Q,V,L}
    es = STI.child_indices_extents(m.node)
    if span_inside(position_span(m.node), m.cells)
        L || return es
        return _shifted(es, first(m.cells) - 1)
    end
    return _members(es, m.cells, Val(L))
end
# Materialise, then rewrite in place: the copy is this call's own, whatever
# container the wrapped leaf hands back.
function _shifted(es, off::Int)
    v = collect(es)
    for i in eachindex(v)
        p, e = @inbounds v[i]
        @inbounds v[i] = (p - off, e)
    end
    return v
end
function _members(es, cells, ::Val{L}) where {L}
    v = collect(es)
    n = 0
    for pe in v
        k = _memberindex(cells, pe[1])
        k == 0 && continue
        n += 1
        @inbounds v[n] = (L ? k : pe[1], pe[2])
    end
    return resize!(v, n)
end

function subtree_indices(m::MaskedCursor{Q,V,L}) where {Q,V,L}
    if span_inside(position_span(m.node), m.cells)
        ps = subtree_indices(m.node)
        L || return ps
        return ps .- (first(m.cells) - 1)
    end
    return _masked_indices!(Int[], m)
end
# Filter the wrapped node's positions where that is cheap against the members
# its span may hold; descend the restriction where the node is far larger (a
# few scattered members under the root of a lattice), never touching extents.
function _masked_indices!(out, m::MaskedCursor{Q,V,L}) where {Q,V,L}
    sp = position_span(m.node)
    if STI.isleaf(m) || sp === nothing || span_inside(sp, m.cells) ||
       split_weight(m.node) <= 4 * _members_within(span_range(sp), m.cells) + 64
        for p in subtree_indices(m.node)
            k = _memberindex(m.cells, p)
            k == 0 || push!(out, L ? k : p)
        end
    else
        for c in STI.getchild(m)
            _masked_indices!(out, c)
        end
    end
    return out
end
_members_within(sp::UnitRange{Int}, cells::AbstractUnitRange) = length(intersect(sp, cells))
_members_within(sp::UnitRange{Int}, cells) =
    max(0, searchsortedlast(cells, last(sp)) - searchsortedfirst(cells, first(sp)) + 1)

# strictly ascending positions, as `Int`s: the shape every `MaskedCursor` reads
_ascending(cells::AbstractUnitRange{<:Integer}) = Int(first(cells)):Int(last(cells))
function _ascending(cells::AbstractVector{<:Integer})
    v = _intvec(cells)
    return _strictly_ascending(v) ? v : sort!(unique(v))
end

# --- the generic entry points ------------------------------------------------

# A `Subset` of this very grid restricts to its positions; any other subset
# through its parent positions.
function subtree(grid, s::Subset; leaf_capacity=Balanced())
    parent(s) === grid || return subtree(grid, parentpositions(s); leaf_capacity)
    n = _capacity(grid, leaf_capacity)
    t = own_tree(grid, n)
    t === nothing && return fallback_subtree(best_manifold(grid), grid, s.cells, n)
    return MaskedCursor{false}(t, s.cells)
end

# A subset's own tree is its parent's, restricted, emitting the subset's
# positions; a parent with no tree of its own leaves the subset its R-tree.
function treeify(s::Subset, n::Integer)
    t = own_tree(s.parent, n)
    t === nothing && return _default_treeify(s, n)
    return MaskedCursor{true}(t, s.cells)
end
