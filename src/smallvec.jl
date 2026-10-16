# ===========================================================================
# A fixed-capacity Int vector with no heap storage, for the small answers a
# provider computes on the spot: neighbour rings, vertex fans.
# ===========================================================================

"""
    SmallInts{N}(values::NTuple{N,Int}, n) <: AbstractVector{Int}
    SmallInts{N}()

The first `n` of `values`, with no heap storage: the shape of an answer that
is at most `N` positions long (neighbour rings, vertex fans). Immutable:
[`push`](@ref GridInterface.push) returns a new one.
"""
struct SmallInts{N} <: AbstractVector{Int}
    v::NTuple{N,Int}
    n::Int
end
SmallInts{N}() where {N} = SmallInts{N}(ntuple(_ -> 0, Val(N)), 0)
Base.size(s::SmallInts) = (s.n,)
Base.IndexStyle(::Type{<:SmallInts}) = IndexLinear()
Base.@propagate_inbounds Base.getindex(s::SmallInts, i::Int) = (@boundscheck checkbounds(s, i); s.v[i])

@noinline _throw_small(::SmallInts{N}) where {N} = throw(ArgumentError("more than $N entries in a SmallInts{$N}"))

"""
    push(s::SmallInts{N}, x) -> SmallInts{N}

`s` with `x` appended; an `ArgumentError` past the capacity.
"""
@inline function push(s::SmallInts{N}, x::Int) where {N}
    s.n < N || _throw_small(s)
    return SmallInts{N}(Base.setindex(s.v, x, s.n + 1), s.n + 1)
end
