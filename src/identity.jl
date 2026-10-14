# ===========================================================================
# Cell identity: the grid owns a dense order; typed ids are a layer on top.
# ===========================================================================

"""
    abstract type AbstractGrid

Optional convenience supertype for grids. The contract is duck-typed, so a type
that already has a supertype elsewhere joins by defining the verbs and
`GridInterface.isgrid(::Type{MyGrid}) = true`.
"""
abstract type AbstractGrid end
Base.broadcastable(g::AbstractGrid) = Ref(g)

"""
    isgrid(x) -> Bool
    isgrid(::Type) -> Bool

Whether `x` claims the GridInterface contract. Define it on the type.
"""
isgrid(::T) where {T} = isgrid(T)
isgrid(::Type) = false
isgrid(::Type{<:AbstractGrid}) = true

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
    iscellindex(x) / iscellindex(::Type) -> Bool

Whether `x` is a typed cell id (not a bare position). Adapters extend it for
foreign id supertypes.
"""
iscellindex(::T) where {T} = iscellindex(T)
iscellindex(::Type) = false
iscellindex(::Type{<:AbstractCellIndex}) = true

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
# ---------------------------------------------------------------------------

"""
    abstract type GridRelation

What "position in the parent" means for a grid. [`NoParent`](@ref) (a complete
grid), [`ParentGrid`](@ref) (a subset: answers `parent` and `globalindex`), or
[`ParentSystem`](@ref) (a multi-level set: answers `system`, not
`globalindex`). Name chosen here; the PRD lists the three relations without
naming a trait.
"""
abstract type GridRelation end
struct NoParent <: GridRelation end
struct ParentGrid <: GridRelation end
struct ParentSystem <: GridRelation end

gridrelation(grid) = NoParent()

"""
    system(grid)

The hierarchical system a grid's cells are drawn from (level grids, subsets of
a level, multi-level sets), or `nothing`.
"""
system(grid) = nothing

"""
    globalindex(grid, c) -> Union{Int,Nothing}

Position of id `c` in the complete grid: for a complete grid the same as
[`localindex`](@ref); for a subset, the position in its parent's complete grid.
A multi-level set refuses (its members have no single parent position).
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
"Local positions `1:ncells(grid)`."
struct Local <: IndexSpace end
"Positions in the complete grid (parent space for a subset)."
struct Global <: IndexSpace end

"""
    Index(Local()), Index(Global()), Index(T)

The "needs" selector a stencil uses to ask for neighbour indices in one space:
local positions, complete-grid positions, or typed ids of scheme `T`.
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
# A subset already stores its parent positions: no id round trip, no `Nothing`.
_toglobal(::NoParent, grid, i) = Int(i)
_toglobal(::ParentGrid, grid, i) = globalindex(grid, cellindex(grid, i))   # provider subsets (DGG PartialGrid)
_toglobal(::ParentSystem, grid, i) = _globalindex(ParentSystem(), grid, i)
toindex(grid, i::Integer, ::Index{Type{T}}) where {T} = cellindex(grid, i, T)
