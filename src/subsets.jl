# ===========================================================================
# Subsets and chunks: a subset is a grid that remembers its parent; a chunking
# is a declared partition into subsets.
# ===========================================================================

"""
    Subset(parent, cells; root = nothing)

The cells of `parent` at the strictly ascending positions `cells` (any
`AbstractVector{Int}`; a `UnitRange` stays a range, a lazy vector stays lazy).
A grid with its own dense order `1:length(cells)`; geometry, identity and traits
forward to the parent. `root` optionally names a hierarchy cell whose
descendants contain every member: it makes `grid_extent` O(1) and lets the tree
start at that cell.

Three index spaces: local (`1:ncells(s)`), parent ([`parentpositions`](@ref)),
complete ([`globalindex`](@ref)).
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
            _strictly_ascending(cells) ||
                throw(ArgumentError("subset positions must be strictly ascending"))
        end
        return new{G,V,R}(parent, cells, root)
    end
end

_strictly_ascending(::AbstractUnitRange) = true
function _strictly_ascending(v::AbstractVector)
    @inbounds for k in 2:length(v)
        v[k] > v[k-1] || return false
    end
    return true
end
function _strictly_ascending(v::RangeCat)
    rs = v.ranges
    for k in 2:length(rs)
        first(rs[k]) > last(rs[k-1]) || return false
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

The positions of a subset's members in its **immediate** parent, in its own
order (the parent space; `subtree` emits these). For a subset of a subset this
differs from [`globalindex`](@ref), which answers in the complete grid.
"""
parentpositions(s::Subset) = s.cells
# provider subsets (a DGG PartialGrid): lazily, through the ids
parentpositions(grid) = _ppositions(gridrelation(grid), grid)
_ppositions(::ParentGrid, grid) = FuncVec(i -> localindex(parent(grid), cellindex(grid, i))::Int, ncells(grid))
_ppositions(r::GridRelation, grid) = throw(ArgumentError("$(typeof(grid)) has no parent grid ($(typeof(r)))"))

gridrelation(::Subset) = ParentGrid()
system(s::Subset) = system(s.parent)
level(s::Subset) = level(s.parent)
GOC.best_manifold(s::Subset) = best_manifold(s.parent)
for f in (:cellshape, :tessellation, :location_is_analytic,
          :neighbors_are_analytic, :cellindextypes, :edgekind, :cell_geometry_is_expensive)
    @eval $f(s::Subset) = $f(s.parent)
end
# one method per core hint as well, so these beat the core's `(Any, Balanced)` forms
for H in (:LeafHint, :Balanced, :ForSweeps, :ForQueries)
    @eval leaf_capacity(s::Subset, h::$H) = leaf_capacity(s.parent, h)
end
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
_ismember(s::Subset, p::Int) = _local_of_parent(s, p) !== nothing

function localindex(s::Subset, c)
    p = localindex(s.parent, c)
    p === nothing && return nothing
    return _local_of_parent(s, p)
end
function localindex(s::Subset, c::Integer)
    p = localindex(s.parent, c)
    p === nothing && return nothing
    return _local_of_parent(s, p)
end
globalindex(s::Subset, c) = localindex(s, c) === nothing ? nothing : globalindex(s.parent, c)

# A subset already stores its parent positions: no id round trip, no `Nothing`.
_toglobal(::ParentGrid, s::Subset, i) = toindex(s.parent, _ppos(s, i), Index(Global()))

for f in (:cell_boundary, :cell_geometry, :cell_centroid, :cell_area, :cell_extent, :cell_corners)
    @eval $f(s::Subset, i::Integer) = $f(s.parent, _ppos(s, i))
end
cell_boundary!(buf::AbstractVector, s::Subset, i::Integer) = cell_boundary!(buf, s.parent, _ppos(s, i))
cell_contains(s::Subset, i::Integer, p) = cell_contains(s.parent, _ppos(s, i), p)

# Tessellating parent: locate in the parent (its tree is shared and cached) and
# map. Overlapping parent: the lowest containing *member*, from the subset's tree.
function cellat(s::Subset, p)
    tessellation(s.parent) isa Overlapping &&
        return fallback_cellat(best_manifold(s), s, check_point(best_manifold(s), p))
    j = cellat(s.parent, p)
    j === nothing && return nothing
    return _local_of_parent(s, j)
end

"""
    grid_extent(s::Subset)

Cheap by contract. With a `root`: `node_extent(system(s), root)`, O(1). On a
subset of a sorted-subtree level grid: the `node_extent` of the lowest common
ancestor of the first and last member, O(depth). Otherwise the union of the
members' `cell_extent`s, O(n) with no tree built.
"""
function grid_extent(s::Subset)
    s.root === nothing || return node_extent(system(s), s.root)
    sys = system(s)
    if sys !== nothing && level(s) !== nothing && has_sorted_subtrees(sys) && ncells(s) > 0
        a = _common_ancestor(sys, s)
        a === nothing || return node_extent(sys, a)
        return _roots_extent(sys, s)
    end
    return _union_extent(s)
end
# members under several roots: the union of the roots they span, O(roots)
function _roots_extent(sys, s)
    ks = _subset_keys(s); lo, hi = first(ks), last(ks); l = level(s)
    e = nothing
    for r in rootcells(sys)
        d = descendant_range(sys, r, l)
        (last(d) < lo || first(d) > hi) && continue
        x = node_extent(sys, r)
        e = e === nothing ? x : Extents.union(e, x)
    end
    return e
end
function _common_ancestor(sys, s)
    a = cellindex(s, 1)
    l = level(s)
    hi = last(_subset_keys(s))
    while level(a) > first(levels(sys)) && last(descendant_range(sys, a, l)) < hi
        a = Base.parent(sys, a)
    end
    last(descendant_range(sys, a, l)) >= hi && return a
    # members under several roots: no single covering cell
    return nothing
end
function _union_extent(g)
    e = cell_extent(g, 1)
    for i in 2:ncells(g)
        e = Extents.union(e, cell_extent(g, i))
    end
    return e
end

# ---------------------------------------------------------------------------
# Subset algebra (same parent), in parent space.
# ---------------------------------------------------------------------------

@noinline _throw_parent() = throw(ArgumentError("subset algebra needs subsets of the same parent"))
_sameparent(a::Subset, b::Subset) = (a.parent === b.parent || a.parent == b.parent) || _throw_parent()

function _merge(op, a::AbstractVector{Int}, b::AbstractVector{Int})
    out = Int[]
    i, j = 1, 1
    na, nb = length(a), length(b)
    @inbounds while i <= na || j <= nb
        if j > nb || (i <= na && a[i] < b[j])
            op === :union || op === :setdiff ? push!(out, a[i]) : nothing
            i += 1
        elseif i > na || b[j] < a[i]
            op === :union && push!(out, b[j])
            j += 1
        else
            op === :union || op === :intersect ? push!(out, a[i]) : nothing
            i += 1; j += 1
        end
    end
    return out
end

"""
    union(a::Subset, b::Subset) / intersect(a, b) / setdiff(a, b) -> Subset

Set algebra on subsets of the same parent, in parent space (a merge of the
sorted parent positions; the result is a new `Subset` of that parent).
"""
Base.union(a::Subset, b::Subset) = (_sameparent(a, b); Subset(a.parent, _merge(:union, a.cells, b.cells)))
Base.intersect(a::Subset, b::Subset) = (_sameparent(a, b); Subset(a.parent, _merge(:intersect, a.cells, b.cells)))
Base.setdiff(a::Subset, b::Subset) = (_sameparent(a, b); Subset(a.parent, _merge(:setdiff, a.cells, b.cells)))

"""
    grow(s::Subset, k = 1; connectivity = Vertex()) -> Subset

`s` plus every parent cell within `k` steps of a member: `union(s, halo(s, k))`.
"""
grow(s::Subset, k::Integer=1; connectivity=Vertex()) =
    Subset(s.parent, _merge(:union, s.cells, halo(s, k; connectivity)))

# ---------------------------------------------------------------------------
# Chunks: a partition into subsets, each with a cheap covering extent.
# ---------------------------------------------------------------------------

"""
    chunks(grid) -> Union{Nothing, AbstractVector of subsets}
    chunks(grid, l) -> AbstractVector of subsets

`chunks(grid)` is the grid's **storage partition** (the tiles or blocks its data
lives in), or `nothing` when it declares none (the default). Consumers coarsen
it by grouping, or ask a hierarchical grid for `chunks(grid, l)`: the level-`l`
blocks of a level grid of a sorted-subtree system, of a subset of one, or of a
multi-level set (a member coarser than `l` is its own chunk). Chunk
collections are lazy (no per-chunk storage on complete level grids). Every
chunk answers [`grid_extent`](@ref) cheaply (stored or O(1)) with an extent
covering its cells, and [`chunkat`](@ref) maps a position to its chunk.

Chunked results are best assembled as local blocks plus a position map
(`parentpositions(chunk)`), not as full-size global matrices.
"""
chunks(grid) = nothing
chunks(grid, l::Integer) = hierarchy_chunks(grid, l)

"""
    chunkat(grid, i) -> Int
    chunkat(chunks, i) -> Int

The index in `chunks(grid)` of the chunk holding position `i`. O(1) or one
`ancestor` for chunks by level; a binary search for window chunks; generic: a
scan with a membership test per chunk.
"""
chunkat(grid, i::Integer) = _chunkat(chunks(grid), grid, Int(i))
_chunkat(::Nothing, grid, i) = throw(ArgumentError("$(typeof(grid)) declares no chunks"))
_chunkat(ch, grid, i) = chunkat(ch, i)
function chunkat(ch::AbstractVector, i::Integer)
    for (k, c) in pairs(ch)
        _ismember(c, Int(i)) && return k
    end
    throw(ArgumentError("position $i is in no chunk"))
end
_ismember(c, p::Int) = p in parentpositions(c)
