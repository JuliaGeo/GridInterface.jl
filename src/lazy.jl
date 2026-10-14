# ===========================================================================
# Lazy, storage-free vectors. Every verb that would otherwise materialise
# O(ncells) ids, positions or chunks returns one of these instead.
# ===========================================================================

"""
    MappedVec(f, v)

Lazy `map(f, v)`: an `AbstractVector` that applies `f` on `getindex`. No storage
beyond `v`.
"""
struct MappedVec{T,F,V<:AbstractVector} <: AbstractVector{T}
    f::F
    v::V
end
MappedVec(f::F, v::V) where {F,V<:AbstractVector} = MappedVec{Base.promote_op(f, eltype(V)),F,V}(f, v)
Base.size(m::MappedVec) = size(m.v)
Base.IndexStyle(::Type{<:MappedVec}) = IndexLinear()
Base.@propagate_inbounds Base.getindex(m::MappedVec, i::Int) = m.f(m.v[i])

"""
    FuncVec{T}(f, n)

Lazy `[f(i) for i in 1:n]`.
"""
struct FuncVec{T,F} <: AbstractVector{T}
    f::F
    n::Int
end
FuncVec(f::F, n::Integer) where {F} = FuncVec{Base.promote_op(f, Int),F}(f, Int(n))
Base.size(m::FuncVec) = (m.n,)
Base.IndexStyle(::Type{<:FuncVec}) = IndexLinear()
Base.@propagate_inbounds function Base.getindex(m::FuncVec, i::Int)
    @boundscheck 1 <= i <= m.n || throw(BoundsError(m, i))
    return m.f(i)
end

"""
    Unsentinel(v, sentinel = 0)

Lazy view of `v` with every `sentinel` entry dropped (the clipped form of a
neighbour list). `length` and `getindex` scan `v`; iteration is one pass. Meant
for short vectors (stencils).
"""
struct Unsentinel{V<:AbstractVector{Int}} <: AbstractVector{Int}
    v::V
    n::Int
end
Unsentinel(v::AbstractVector{Int}) = Unsentinel(v, count(!=(0), v))
Base.size(u::Unsentinel) = (u.n,)
Base.IndexStyle(::Type{<:Unsentinel}) = IndexLinear()
function Base.getindex(u::Unsentinel, i::Int)
    @boundscheck 1 <= i <= u.n || throw(BoundsError(u, i))
    k = 0
    for x in u.v
        x == 0 && continue
        k += 1
        k == i && return x
    end
    error("unreachable")
end
function Base.iterate(u::Unsentinel, s::Int=1)
    v = u.v
    while s <= length(v)
        @inbounds x = v[s]
        x != 0 && return (x, s + 1)
        s += 1
    end
    return nothing
end

"""
    RangeCat(ranges, offsets)

Lazy concatenation of ascending, disjoint `UnitRange`s: the positions of an
expanded multi-level set at its reference level. `getindex` is a binary search
over the cumulative counts, so storage is O(number of ranges), not O(cells).
"""
struct RangeCat{R<:AbstractVector{UnitRange{Int}}} <: AbstractVector{Int}
    ranges::R
    offsets::Vector{Int}      # offsets[k] = cells before range k; offsets[end] = total
end
function RangeCat(ranges::AbstractVector{UnitRange{Int}})
    offs = Vector{Int}(undef, length(ranges) + 1)
    offs[1] = 0
    for (k, r) in enumerate(ranges)
        offs[k+1] = offs[k] + length(r)
    end
    return RangeCat(ranges, offs)
end
Base.size(r::RangeCat) = (r.offsets[end],)
Base.IndexStyle(::Type{<:RangeCat}) = IndexLinear()
function Base.getindex(r::RangeCat, i::Int)
    @boundscheck 1 <= i <= r.offsets[end] || throw(BoundsError(r, i))
    k = searchsortedlast(r.offsets, i - 1)
    return first(r.ranges[k]) + (i - 1 - r.offsets[k])
end
Base.issorted(::RangeCat) = true
