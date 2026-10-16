# ===========================================================================
# Cell identity: the grid owns a dense order; typed ids are a layer on top.
# ===========================================================================

"""
    abstract type AbstractGrid

Optional convenience supertype for grids. The contract is duck-typed, so a type
that already has a supertype elsewhere joins by defining the verbs and
`GridInterface.isgrid(::MyGrid) = true`.
"""
abstract type AbstractGrid end
Base.broadcastable(g::AbstractGrid) = Ref(g)

"""
    isgrid(x) -> Bool

Whether `x` claims the GridInterface contract. An instance trait (never
`Type{<:X}`): `true` for an [`AbstractGrid`](@ref) or a wrapper, `false` by
default; a foreign grid type adds `isgrid(::MyGrid) = true`.
"""
isgrid(x) = false
isgrid(::AbstractGrid) = true

"""
    abstract type AbstractCellIndex

Supertype of typed cell ids. Contract: isbits, a canonical total order
(`isless`, consistent `==`/`hash`), and, for an id that belongs to a
hierarchical system, [`level(c)`](@ref level) (a binding law: the level is how a
mixed-level collection finds the grid that owns a cell's geometry).

Ids of existing packages that cannot subtype this (DGG's ids subtype
`DiscreteGlobalGrids.AbstractCellIndex`) must add the typed-id geometry forms
(`cell_boundary(grid, c)`, ...) themselves: the generic typed-id forms dispatch
on this supertype, not on `iscellindex`.
"""
abstract type AbstractCellIndex end

"""
    iscellindex(c) -> Bool

Whether `c` is a typed cell id (not a bare position): an
[`AbstractCellIndex`](@ref) by default. Adapters add methods for foreign id
supertypes (`iscellindex(::ForeignId) = true`).
"""
iscellindex(c) = c isa AbstractCellIndex

"""
    ncells(grid) -> Int

Length of the grid's dense order. **Required.** O(1), constant for the
lifetime of the grid object.
"""
function ncells end

"""
    level(c) -> Int
    level(grid) -> Union{Int,Nothing}

The level of a hierarchical id (required on every id of a hierarchical system),
or the level a single-level grid is drawn from (`nothing` for standalone grids
and multi-level sets).
"""
function level end
level(grid) = nothing

@noinline _throw_bounds(grid, i) = throw(BoundsError(grid, i))

"""
    cellindex(grid, i::Integer)

The canonical id of the cell at position `i`. A grid with no id scheme names
cells by position, so the default returns `i` after a bounds check.
Out-of-range `i` throws `BoundsError`.
"""
function cellindex(grid, i::Integer)
    1 <= i <= ncells(grid) || _throw_bounds(grid, i)
    return Int(i)
end

"""
    cellindextypes(grid) -> Tuple of types

The id schemes `grid` names cells in, canonical first. Default: the type
`cellindex(grid, 1)` returns (`Int` for an id-less grid).
"""
cellindextypes(grid) = ncells(grid) > 0 ? (typeof(cellindex(grid, 1)),) : (Int,)

@noinline _throw_scheme(grid, T) = throw(ArgumentError(
    "id scheme $T is not declared by $(typeof(grid)); cellindextypes = $(cellindextypes(grid))"))

_declares(grid, ::Type{T}) where {T} = any(S -> S === T, cellindextypes(grid))

"""
    cellindex(grid, i::Integer, T::Type)

The id of the cell at position `i` in scheme `T`. An undeclared `T` is refused
with an `ArgumentError` before any work.
"""
function cellindex(grid, i::Integer, ::Type{T}) where {T}
    _declares(grid, T) || _throw_scheme(grid, T)
    c = cellindex(grid, i)
    c isa T && return c
    return reindex(T, grid, c)
end

"""
    reindex(T, grid, c)

Convert id `c` to scheme `T` (same cell). `reindex(typeof(c), grid, c) === c`.
Undeclared `T` throws `ArgumentError`; a grid that declares several schemes
implements the conversion.
"""
function reindex(::Type{T}, grid, c) where {T}
    c isa T && return c
    _declares(grid, T) || _throw_scheme(grid, T)
    throw(ArgumentError("$(typeof(grid)) declares scheme $T but implements no `reindex` to it"))
end

"""
    localindex(grid, c) -> Union{Int,Nothing}

Position of the cell with id `c` (any declared scheme; canonicalised first), or
`nothing` for a non-member. For an id-less grid this is a bounds check. The
generic method for typed ids is a linear scan; every grid with an id scheme
should override it (the conformance suite measures it).
"""
function localindex(grid, c::Integer)
    T = first(cellindextypes(grid))
    T <: Integer || return _scan_localindex(grid, c)
    return 1 <= c <= ncells(grid) ? Int(c) : nothing
end
function localindex(grid, c)
    T = first(cellindextypes(grid))
    if !(c isa T)
        _declares(grid, typeof(c)) || return nothing
        c = reindex(T, grid, c)
    end
    return _scan_localindex(grid, c)
end
function _scan_localindex(grid, c)
    for i in 1:ncells(grid)
        cellindex(grid, i) == c && return i
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Relations a collection can have: none, a parent grid, a parent system.
# Three index spaces follow from them:
#   local     1:ncells(grid)                    every grid
#   parent    parentpositions(s)[i]             a subset, in the grid it was cut from
#   complete  globalindex(s, c)                 a subset, in the complete grid at the
#                                               root of its parent chain
# For a subset of a complete grid, parent and complete coincide. For a subset of
# a subset they differ: chunk assembly into the subset's own matrix uses the
# parent space, halo exchange with the full grid uses the complete space.
# ---------------------------------------------------------------------------

"""
    abstract type GridRelation

What "position in the parent" means for a grid: [`NoParent`](@ref),
[`ParentGrid`](@ref) or [`ParentSystem`](@ref), answered by
[`gridrelation`](@ref).
"""
abstract type GridRelation end
"A complete grid: local positions are the only positions."
struct NoParent <: GridRelation end
"A subset of a parent grid: answers `parent`, [`parentpositions`](@ref) and [`globalindex`](@ref)."
struct ParentGrid <: GridRelation end
"A multi-level set: answers [`system`](@ref), refuses [`globalindex`](@ref)."
struct ParentSystem <: GridRelation end

"""
    gridrelation(grid) -> GridRelation

`NoParent()` (default), `ParentGrid()` or `ParentSystem()`.
"""
gridrelation(grid) = NoParent()

"""
    system(grid)

The hierarchical system a grid's cells are drawn from, or `nothing`.
**Obligation** for the level grids of a hierarchical system (together with
[`level(grid)`](@ref level)): it is how the generic hierarchy cursor, chunks
and multi-level sets find the hierarchy. Subsets and multi-level sets answer it
too.
"""
system(grid) = nothing

"""
    globalindex(grid, c) -> Union{Int,Nothing}

Position of id `c` in the **complete** grid: for a complete grid the same as
[`localindex`](@ref); for a subset (at any nesting depth) the position in the
complete grid at the root of its parent chain. For the position in the
immediate parent use [`parentpositions`](@ref). A multi-level set refuses (its
members have no single complete-grid position).
"""
globalindex(grid, c) = _globalindex(gridrelation(grid), grid, c)
_globalindex(::NoParent, grid, c) = localindex(grid, c)
_globalindex(::ParentGrid, grid, c) =
    localindex(grid, c) === nothing ? nothing : globalindex(parent(grid), c)
_globalindex(::ParentSystem, grid, c) = throw(ArgumentError(
    "$(typeof(grid)) has a parent system, not a parent grid: `globalindex` has no meaning. " *
    "Use the typed id (`cellindex`) or `descendant_range` at a reference level."))

# ---------------------------------------------------------------------------
# Stencil selectors: which index space a consumer wants neighbours in.
# ---------------------------------------------------------------------------

abstract type IndexSpace end
"""
    Local()

Index space selector: local positions `1:ncells(grid)`. On a subset,
`neighbors(s, i, Index(Local()))` keeps every parent neighbour slot and puts
`0` (the outside sentinel, [`OUTSIDE`](@ref)) where the neighbour is not a member.
"""
struct Local <: IndexSpace end
"""
    Global()

Index space selector: positions in the complete grid. On a subset,
`neighbors(s, i, Index(Global()))` is the **unclipped** neighbourhood: every
neighbour, inside or outside the subset, in complete-grid positions.
"""
struct Global <: IndexSpace end

"""
    OUTSIDE = 0

The sentinel `neighbors(s, i, Index(Local()))` puts in a slot whose neighbour
lies outside subset `s`.
"""
const OUTSIDE = 0

"""
    Index(Local()), Index(Global()), Index(T)

The "needs" selector a stencil uses to ask for neighbour indices in one space:
local positions (with [`OUTSIDE`](@ref) for non-members), complete-grid
positions, or typed ids of scheme `T`. The selector forms never clip: slot `k`
is the same neighbour in every space.
"""
struct Index{S}
    space::S
end
Index(::Type{T}) where {T} = Index{Type{T}}(T)

"""
    toindex(grid, i::Integer, sel::Index)

Position `i` of `grid` expressed in the space `sel` selects.
"""
toindex(grid, i::Integer, ::Index{Local}) = Int(i)
toindex(grid, i::Integer, ::Index{Global}) = _toglobal(gridrelation(grid), grid, i)::Int
_toglobal(::NoParent, grid, i) = Int(i)
_toglobal(::ParentGrid, grid, i) = globalindex(grid, cellindex(grid, i))   # provider subsets (DGG PartialGrid)
_toglobal(::ParentSystem, grid, i) = _globalindex(ParentSystem(), grid, i)
toindex(grid, i::Integer, ::Index{Type{T}}) where {T} = cellindex(grid, i, T)
