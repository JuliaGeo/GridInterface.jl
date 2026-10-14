# ===========================================================================
# Hierarchical grid systems: the primitives. The derived analysis layer stays
# in the provider (DGG).
# ===========================================================================

"""
    abstract type AbstractHierarchicalGridSystem

A family of grids, one per level, with parent/child relations between cells of
adjacent levels. Optional supertype; foreign systems join with
`GridInterface.issystem(::Type{S}) = true`.

Primitives (all on typed ids `c`, every one of which answers `level(c)`):

| verb | contract |
|---|---|
| `levels(sys)` | the level numbers (an ascending range) |
| `rootcells(sys)` | ids at the coarsest level, ascending |
| `parent(sys, c)` | `Base.parent`; the analytic parent |
| `children(sys, c)` | ascending canonical order |
| `levelgrid(sys, l)` | a complete grid for level `l`; **one concrete type for every level** |
| `has_sorted_subtrees(sys)` | `true` obliges an exact `descendant_range` |
| `descendant_range(sys, c, l)` | positions in `levelgrid(sys, l)` of `c`'s descendants |
| `node_extent(sys, c)` | covering law: contains the geometry of every descendant |
| `has_congruent_refinement(sys)` | children tile their parent exactly |
"""
abstract type AbstractHierarchicalGridSystem end
Base.broadcastable(s::AbstractHierarchicalGridSystem) = Ref(s)

issystem(::T) where {T} = issystem(T)
issystem(::Type) = false
issystem(::Type{<:AbstractHierarchicalGridSystem}) = true

function levels end
function rootcells end
function children end
function levelgrid end
function node_extent end

has_sorted_subtrees(sys) = false
has_congruent_refinement(sys) = false

"""
    descendant_range(sys, c, l) -> UnitRange{Int}

Positions in `levelgrid(sys, l)` of the level-`l` descendants of `c`
(`l >= level(c)`). Obliged and exact when `has_sorted_subtrees(sys)`.
"""
function descendant_range end

"""
    ancestor(sys, c, l)

The level-`l` ancestor of `c` (`l <= level(c)`). Generic: the parent loop.
"""
function ancestor(sys, c, l::Integer)
    lc = level(c)
    l <= lc || throw(ArgumentError("ancestor level $l is below the cell's level $lc"))
    while lc > l
        c = Base.parent(sys, c)
        lc -= 1
    end
    return c
end

"""
    hierarchy_chunks(grid, l) -> Vector{Subset}

The partition of a complete level grid of a sorted-subtree system into the
descendant blocks of every level-`l` cell: each chunk is a `Subset` over a
`UnitRange` rooted at its ancestor, so its extent is `node_extent(sys, root)`.
"""
function hierarchy_chunks(grid, l::Integer)
    sys = system(grid)
    lg = level(grid)
    (sys === nothing || lg === nothing) && throw(ArgumentError("not a level grid"))
    has_sorted_subtrees(sys) || throw(ArgumentError("hierarchy_chunks needs sorted subtrees"))
    l <= lg || throw(ArgumentError("chunk level $l is finer than the grid's level $lg"))
    out = Subset{typeof(grid),UnitRange{Int},eltype(_cells_at(sys, l))}[]
    # A subset of a level grid (a `Subset`, a DGG `PartialGrid`) numbers its own
    # cells: map each descendant range of the parent level to the window of the
    # subset's ascending parent positions it covers, and drop empty windows.
    complete = gridrelation(grid) isa NoParent
    ks = complete ? nothing : _subset_keys(grid)
    for c in _cells_at(sys, l)
        r = descendant_range(sys, c, lg)
        if !complete
            r = searchsortedfirst(ks, first(r)):searchsortedlast(ks, last(r))
            isempty(r) && continue
        end
        push!(out, Subset(grid, r; root=c))
    end
    return out
end

function _cells_at(sys, l)
    g = levelgrid(sys, l)
    return [cellindex(g, i) for i in 1:ncells(g)]
end
