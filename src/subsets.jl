# ===========================================================================
# Subsets and chunks: a subset is a grid that remembers its parent; a chunking
# is a declared partition into subsets.
# ===========================================================================

"""
    Subset(parent, cells; root = nothing)

The cells of `parent` at the strictly ascending positions `cells` (any
`AbstractVector{Int}`; a `UnitRange` stays a range). A grid with its own dense
order `1:length(cells)`; geometry, identity and traits forward to the parent.
`root` optionally names a hierarchy cell whose descendants at the parent's
level contain every member: it makes `grid_extent` exact and lets the tree start
at that cell.
"""
struct Subset{G,V<:AbstractVector{Int},R} <: AbstractGrid
    parent::G
    cells::V
    root::R
    function Subset(parent::G, cells::V; root::R=nothing) where {G,V<:AbstractVector{Int},R}
        n = ncells(parent)
        if !isempty(cells)
            (1 <= first(cells) && last(cells) <= n) ||
                throw(ArgumentError("subset positions must lie in 1:$n"))
            cells isa AbstractUnitRange || _strictly_ascending(cells) ||
                throw(ArgumentError("subset positions must be strictly ascending"))
        end
        return new{G,V,R}(parent, cells, root)
    end
end

function _strictly_ascending(v)
    @inbounds for k in 2:length(v)
        v[k] > v[k-1] || return false
    end
    return true
end

"""
    subset(grid, cells)

The subset of `grid` at positions `cells` (strictly ascending positions or a
range). The representation is free: providers return their own subset types
(a raster tile, a DGG `PartialGrid`); the generic one is [`Subset`](@ref).
Typed ids are accepted and resolved through `localindex`.
"""
subset(grid, cells::AbstractVector{<:Integer}) = Subset(grid, _intvec(cells))
subset(grid, cells::AbstractVector) =
    Subset(grid, sort!([_pos(grid, c) for c in cells]))
_intvec(v::AbstractVector{Int}) = v
_intvec(v::AbstractVector{<:Integer}) = Int.(v)

Base.parent(s::Subset) = s.parent

"""
    parentpositions(s) -> AbstractVector{Int}

The parent positions of a subset's members, in its own order.
"""
parentpositions(s::Subset) = s.cells

gridrelation(::Subset) = ParentGrid()
system(s::Subset) = system(s.parent)
level(s::Subset) = level(s.parent)
GOC.best_manifold(s::Subset) = best_manifold(s.parent)
cellshape(s::Subset) = cellshape(s.parent)
tessellation(s::Subset) = tessellation(s.parent)
extent_is_expensive(s::Subset) = extent_is_expensive(s.parent)
location_is_analytic(s::Subset) = location_is_analytic(s.parent)
neighbors_are_analytic(s::Subset) = neighbors_are_analytic(s.parent)
cellindextypes(s::Subset) = cellindextypes(s.parent)
reindex(::Type{T}, s::Subset, c) where {T} = reindex(T, s.parent, c)

ncells(s::Subset) = length(s.cells)

Base.@propagate_inbounds function _ppos(s::Subset, i::Integer)
    @boundscheck (1 <= i <= length(s.cells)) || _throw_bounds(s, i)
    return @inbounds s.cells[i]
end

cellindex(s::Subset, i::Integer) = cellindex(s.parent, _ppos(s, i))

# parent position -> own position
@inline function _local_of_parent(s::Subset, p::Int)
    cs = s.cells
    if cs isa AbstractUnitRange
        return p in cs ? p - first(cs) + 1 : nothing
    end
    k = searchsortedfirst(cs, p)
    return (k <= length(cs) && @inbounds(cs[k]) == p) ? k : nothing
end

function localindex(s::Subset, c)
    p = localindex(s.parent, c)
    p === nothing && return nothing
    return _local_of_parent(s, p)
end
localindex(s::Subset, c::Integer) = _subset_localindex_int(s, c)
function _subset_localindex_int(s, c)
    p = localindex(s.parent, c)
    p === nothing && return nothing
    return _local_of_parent(s, p)
end

for f in (:cell_boundary, :cell_geometry, :cell_centroid, :cell_area, :cell_extent)
    @eval $f(s::Subset, i::Integer) = $f(s.parent, _ppos(s, i))
end
cell_contains(s::Subset, i::Integer, p) = cell_contains(s.parent, _ppos(s, i), p)

function cellat(s::Subset, p)
    j = cellat(s.parent, p)
    j === nothing && return nothing
    return _local_of_parent(s, j)
end

function grid_extent(s::Subset)
    s.root === nothing && return STI.node_extent(treeify(s))
    return node_extent(system(s), s.root)
end

# ---------------------------------------------------------------------------
# Chunks: a partition into subsets, each with an extent. Nothing more.
# ---------------------------------------------------------------------------

"""
    chunks(grid) -> Union{Nothing, AbstractVector of subsets}

A declared partition of `grid` into subsets, each answering
[`grid_extent`](@ref). `nothing` when the grid declares none (the default).
Level grids of systems with sorted subtrees can use
[`hierarchy_chunks`](@ref).
"""
chunks(grid) = nothing

# A subset already stores its parent positions: no id round trip, no `Nothing`.
_toglobal(::ParentGrid, s::Subset, i) = toindex(s.parent, s.cells[i], Index(Global()))
