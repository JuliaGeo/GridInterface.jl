# ===========================================================================
# Hierarchical grid systems: the primitives. The derived analysis layer stays
# in the provider (DGG).
# ===========================================================================

"""
    abstract type AbstractHierarchicalGridSystem

A family of grids, one per level, with parent/child relations between cells of
adjacent levels. Optional supertype; foreign systems join with
`GridInterface.issystem(::S) = true`.

Obligations (all on typed ids `c`; the ids subtype [`AbstractCellIndex`](@ref)):

| verb | contract |
|---|---|
| `level(c)` | the id's level |
| `levels(sys)` | the level numbers (an ascending range) |
| `rootcells(sys)` | ids at the coarsest level, ascending |
| `parent(sys, c)` | `Base.parent`; the analytic parent |
| `children(sys, c)` | ascending canonical order |
| `levelgrid(sys, l)` | a complete grid for level `l`; **one concrete type for every level** |
| `system(levelgrid(sys, l))`, `level(levelgrid(sys, l))` | `sys` and `l`: how the generic cursor, chunks and multi-level sets find the hierarchy |
| `localindex(levelgrid(sys, l), c)` | analytic (no search) |
| `isless(a, b)` | matches position order within a level |
| `has_sorted_subtrees(sys)` | `true` obliges an exact `descendant_range` |
| `descendant_range(sys, c, l)` | positions in `levelgrid(sys, l)` of `c`'s descendants |
| `node_extent(sys, c)` | covering law: contains the geometry of every descendant |
| `Extents.union(e1, e2)` | defined on the `node_extent` type |
| `has_congruent_refinement(sys)` | children tile their parent exactly |

Radius rule: every level grid's `best_manifold` carries the same radius, which
is the unit [`cell_area`](@ref) is measured in.
"""
abstract type AbstractHierarchicalGridSystem end
Base.broadcastable(s::AbstractHierarchicalGridSystem) = Ref(s)

"""
    issystem(x) -> Bool

Whether `x` claims the hierarchical-system contract. An instance trait: `true`
for an [`AbstractHierarchicalGridSystem`](@ref), `false` by default; a foreign
system adds `issystem(::S) = true`.
"""
issystem(x) = false
issystem(::AbstractHierarchicalGridSystem) = true

"""
    levels(sys) -> AbstractUnitRange{Int}

The system's level numbers, ascending.
"""
function levels end

"""
    rootcells(sys)

The ids at the coarsest level, ascending (an indexable collection).
"""
function rootcells end

"""
    children(sys, c)

The ids one level below `c` whose parent is `c`, ascending.
"""
function children end

"""
    levelgrid(sys, l)

The complete grid of level `l`. One concrete return type for every `l`, so a
multi-level set forwards geometry type-stably. It answers `system` and `level`.
"""
function levelgrid end

"""
    node_extent(sys, c)

An extent (box or `SphericalCap`) containing the geometry of every descendant
of `c` (the covering law; it need not contain a child's `node_extent`).
"""
function node_extent end

"""
    has_sorted_subtrees(sys) -> Bool

Whether every cell's descendants at any finer level occupy one contiguous
position range ([`descendant_range`](@ref)). Enables the window-mode cursor,
chunks by level, and O(1) `subtree`. Default `false`.
"""
has_sorted_subtrees(sys) = false

"""
    has_congruent_refinement(sys) -> Bool

Whether children tile their parent exactly (then a multi-level set is
tessellating). Default `false`.
"""
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

# ---------------------------------------------------------------------------
# Chunks by level: lazy collections of Subsets, never materialised.
# ---------------------------------------------------------------------------

"""
    LevelChunks(l)

A chunking by level, as a value: `mockchunks(grid, LevelChunks(l))` is `grid`
with `chunks` answering `chunks(grid, l)`, the level-`l` blocks (see
[`chunks`](@ref), [`mockchunks`](@ref)).
"""
struct LevelChunks
    level::Int
end

"""
    AncestorChunks(grid, sys, chunkgrid, lg, n)

What `chunks(grid, l)` returns on a complete level grid of a sorted-subtree
system: chunk `k` is the `Subset` over `descendant_range` of the `k`-th
level-`l` cell, rooted at that cell. A lazy `AbstractVector` with no storage;
`getindex` is O(1) plus one `cellindex`; [`chunkat`](@ref) is one `ancestor`.
"""
struct AncestorChunks{G,S,ID,LG} <: AbstractVector{Subset{G,UnitRange{Int},ID}}
    grid::G
    sys::S
    chunkgrid::LG        # levelgrid(sys, l)
    lg::Int              # the grid's own level
    n::Int
end
Base.size(c::AncestorChunks) = (c.n,)
Base.IndexStyle(::Type{<:AncestorChunks}) = IndexLinear()
function Base.getindex(c::AncestorChunks, k::Int)
    @boundscheck 1 <= k <= c.n || throw(BoundsError(c, k))
    a = cellindex(c.chunkgrid, k)
    return Subset(c.grid, descendant_range(c.sys, a, c.lg); root=a)
end
function chunkat(c::AncestorChunks, i::Integer)
    a = ancestor(c.sys, cellindex(c.grid, i), level(c.chunkgrid))
    return localindex(c.chunkgrid, a)::Int
end

"""
    WindowChunks(grid, starts, roots)

Chunks of a non-complete grid (a subset of a level grid, a multi-level set): chunk
`k` is the `Subset` over own positions `starts[k]:starts[k+1]-1`, rooted at
`roots[k]` (its covering cell). Storage is O(number of chunks).
"""
struct WindowChunks{G,ID} <: AbstractVector{Subset{G,UnitRange{Int},ID}}
    grid::G
    starts::Vector{Int}   # length nchunks + 1
    roots::Vector{ID}
end
Base.size(c::WindowChunks) = (length(c.roots),)
Base.IndexStyle(::Type{<:WindowChunks}) = IndexLinear()
function Base.getindex(c::WindowChunks, k::Int)
    @boundscheck 1 <= k <= length(c.roots) || throw(BoundsError(c, k))
    return Subset(c.grid, c.starts[k]:c.starts[k+1]-1; root=c.roots[k])
end
chunkat(c::WindowChunks, i::Integer) = searchsortedlast(c.starts, Int(i))

"""
    hierarchy_chunks(grid, l)

The level-`l` chunks of `grid` (see [`chunks`](@ref)`(grid, l)`).
"""
function hierarchy_chunks(grid, l::Integer)
    sys = system(grid)
    lg = level(grid)
    (sys === nothing || lg === nothing) && throw(ArgumentError("not a level grid (or a subset of one)"))
    has_sorted_subtrees(sys) || throw(ArgumentError("chunks by level need sorted subtrees"))
    l <= lg || throw(ArgumentError("chunk level $l is finer than the grid's level $lg"))
    if gridrelation(grid) isa NoParent
        cg = levelgrid(sys, l)
        ID = typeof(cellindex(cg, 1))
        return AncestorChunks{typeof(grid),typeof(sys),ID,typeof(cg)}(grid, sys, cg, lg, ncells(cg))
    end
    # A subset of a level grid numbers its own cells: one window per ancestor,
    # found by a binary search over the subset's ascending parent positions.
    ks = _subset_keys(grid)
    n = length(ks)
    ID = typeof(cellindex(grid, 1))
    starts = Int[]; roots = ID[]
    k = 1
    while k <= n
        a = ancestor(sys, cellindex(grid, k), l)
        r = descendant_range(sys, a, lg)
        hi = searchsortedlast(ks, last(r), k, n, Base.Order.Forward)
        push!(starts, k); push!(roots, a)
        k = hi + 1
    end
    push!(starts, n + 1)
    return WindowChunks{typeof(grid),ID}(grid, starts, roots)
end
