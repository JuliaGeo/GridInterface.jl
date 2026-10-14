# DGGAdapter: DiscreteGlobalGrids level grids, partial grids, multi-order sets
# and hierarchical systems satisfy GridInterface. Every method here is defined
# from outside on DGG's own types; DGG is not edited.
#
#   include("prototype/adapters/DGGAdapter.jl"); using .DGGAdapter
#
# What maps to what:
#   DGG.AbstractGrid (HierarchicalLevelGrid, PartialGrid, AuthalicGrid, ...)  → a GridInterface grid
#   DGG.AbstractHierarchicalGridSystem                                        → a GridInterface system
#   DGG.AbstractCellIndex                                                     → typed ids (level(c) binding)
#   DGG.MultiOrderCellSet                                                     → GridInterface.MultiLevelSet (O(1) wrap)
# treeify of a level grid / partial grid is GridInterface's generic hierarchy
# cursor; DGG's own cursor stays reachable as `dgg_tree(g)` (wrap with GridAndTree).
module DGGAdapter

import GridInterface as GrI
using GridInterface
import DiscreteGlobalGrids as DGG
const FB = DGG.Fallbacks
import GeometryOpsCore as GOC

export dgg_multilevel, dgg_tree, dgg_systems, DGG

const G = DGG.AbstractGrid
const SYS = DGG.AbstractHierarchicalGridSystem
const ID = DGG.AbstractCellIndex
const HLG = DGG.HierarchicalLevelGrid

# --- recognition --------------------------------------------------------------
GrI.isgrid(::Type{<:G}) = true
GrI.issystem(::Type{<:SYS}) = true
GrI.iscellindex(::Type{<:ID}) = true
GrI.level(c::ID) = DGG.level(c)

# --- identity -----------------------------------------------------------------
GrI.ncells(g::G) = DGG.ncells(g)
GrI.cellindex(g::G, i::Integer) = DGG.cellindex(g, Int(i))
GrI.cellindex(g::G, i::Integer, ::Type{T}) where {T} =
    (GrI._declares(g, T) || GrI._throw_scheme(g, T); DGG.cellindex(g, Int(i), T))
GrI.cellindextypes(g::G) = DGG.cellindextypes(g)
GrI.reindex(::Type{T}, g::G, c) where {T} = DGG.reindex(T, DGG.system(g), c)
GrI.localindex(g::G, c::ID) = DGG.localindex(g, c)
GrI.localindex(g::G, c::Integer) = nothing            # DGG grids have no position-named ids
GrI.globalindex(g::G, c::ID) = DGG.globalindex(g, c)
GrI.system(g::G) = DGG.system(g)
GrI.level(g::G) = DGG.level(g)

GrI.gridrelation(::DGG.PartialGrid) = GrI.ParentGrid()
Base.parent(pg::DGG.PartialGrid) = pg.complete          # (DGG defines no `parent` for it)

# --- geometry -----------------------------------------------------------------
# Cell-shape class, declared once per system (PRD: per system, holds for every
# level grid). ConvexCells only where the boundary is the great-circle polygon of
# the cell; systems that densify curved chart edges declare SimpleCells.
dgg_cellshape(::Any) = GrI.SimpleCells()
dgg_cellshape(::DGG.S2System) = GrI.ConvexCells()
dgg_cellshape(::DGG.H3System) = GrI.ConvexCells()
dgg_cellshape(::DGG.IGeo7System) = GrI.ConvexCells()
dgg_cellshape(s::DGG.AuthalicSystem) = dgg_cellshape(s.system)
GrI.cellshape(g::G) = dgg_cellshape(DGG.system(g))
GrI.extent_is_expensive(g::G) = !(FB.cell_cap_is_cheap(g) isa Val{true})

@inline _id(g, i) = DGG.cellindex(g, Int(i))
GrI.cell_boundary(g::G, i::Integer) = DGG.cell_boundary(g, _id(g, i))
GrI.cell_geometry(g::G, i::Integer) = DGG.cell_polygon(g, _id(g, i))
GrI.cell_corners(g::G, i::Integer) = DGG.cell_corners(g, _id(g, i))
GrI.cell_centroid(g::G, i::Integer) = DGG.cell_centroid(g, _id(g, i))
GrI.cell_area(g::G, i::Integer) = DGG.cell_area(g, _id(g, i))
GrI.cell_extent(g::G, i::Integer) = FB.cell_cap(g, _id(g, i))
# typed-id forms: no localindex round trip (what a MultiLevelSet forwards through)
GrI.cell_boundary(g::G, c::ID) = DGG.cell_boundary(g, c)
GrI.cell_geometry(g::G, c::ID) = DGG.cell_polygon(g, c)
GrI.cell_centroid(g::G, c::ID) = DGG.cell_centroid(g, c)
GrI.cell_area(g::G, c::ID) = DGG.cell_area(g, c)
GrI.cell_extent(g::G, c::ID) = FB.cell_cap(g, c)

# --- location -----------------------------------------------------------------
function GrI.location_is_analytic(g::G)
    s = DGG.system(g)
    return s !== nothing && DGG.has_direct_location(s)
end
function GrI.cellat(g::G, p)
    GrI.location_is_analytic(g) || return GrI.fallback_cellat(GOC.best_manifold(g), g, p)
    c = DGG.cellat(g, p)
    c === nothing && return nothing
    return DGG.localindex(g, c)
end

# --- neighbours ---------------------------------------------------------------
_dconn(::GrI.Vertex) = DGG.Vertex()
_dconn(::GrI.Edge) = DGG.Edge()
_gwind(::DGG.CounterClockwise) = GrI.CounterClockwise()
_gwind(::DGG.Clockwise) = GrI.Clockwise()
_gwind(::DGG.CustomOrder) = GrI.CustomOrder()
_gwind(::DGG.Unordered) = GrI.Unordered()
GrI.neighbors_are_analytic(g::G) = DGG.system(g) !== nothing
GrI.neighbors(g::G, i::Integer; connectivity::GrI.Connectivity=GrI.Vertex()) =
    DGG.neighbors(g, Int(i), 1; connectivity=_dconn(connectivity))
GrI.ring(g::G, i::Integer, k::Integer; connectivity::GrI.Connectivity=GrI.Vertex()) =
    DGG.ring(g, Int(i), k; connectivity=_dconn(connectivity))
GrI.maxneighbors(g::G, c::GrI.Connectivity) = DGG.maxneighbors(g, _dconn(c))
GrI.winding(g::G, c::GrI.Connectivity) = _gwind(DGG.winding(g, _dconn(c)))

# --- tree and chunks ------------------------------------------------------------
# The generic hierarchy cursor (window mode on sorted systems, selection mode on A5).
GrI.treeify(g::HLG) = GrI.hierarchy_tree(g)
"DGG's own tree for a grid (HierarchicalGridCursor, MemoBlockCursor, ...): `GridAndTree(g, dgg_tree(g))`."
dgg_tree(g) = DGG.treeify(g)

function GrI.chunks(g::HLG)
    s = DGG.system(g)
    DGG.has_sorted_subtrees(s) || return nothing
    l0 = first(DGG.levels(s))
    return GrI.hierarchy_chunks(g, max(l0, g.level - 3))
end

# --- hierarchy primitives --------------------------------------------------------
GrI.levels(s::SYS) = DGG.levels(s)
GrI.rootcells(s::SYS) = DGG.rootcells(s)
GrI.children(s::SYS, c::ID) = DGG.children(s, c)
GrI.levelgrid(s::SYS, l::Integer) = DGG.levelgrid(s, l)
GrI.has_sorted_subtrees(s::SYS) = DGG.has_sorted_subtrees(s)
GrI.has_congruent_refinement(s::SYS) = DGG.has_congruent_refinement(s)
GrI.descendant_range(s::SYS, c::ID, l::Integer) = DGG.descendant_range(s, c, l)
GrI.node_extent(s::SYS, c::ID) = DGG.node_extent(s, c)
GrI.ancestor(s::SYS, c::ID, l::Integer) = DGG.ancestor(s, c, l)
# (`Base.parent(sys, c)` is DGG's own method on Base.parent: nothing to add.)

# --- multi-order sets -------------------------------------------------------------
"""
    dgg_multilevel(set::DGG.MultiOrderCellSet) -> GridInterface.MultiLevelSet

O(1) on sorted systems: DGG's storage (curve order at the reference level, one
key per cell) is exactly the MultiLevelSet's sorted layout.
"""
function dgg_multilevel(set::DGG.MultiOrderCellSet)
    s = set.system
    if DGG.has_sorted_subtrees(s)
        return GrI.MultiLevelSet{typeof(s),eltype(set.cells)}(s, set.cells, set.keys, Int[], set.reference_level, true)
    end
    return GrI.MultiLevelSet(s, set.cells)
end

"The registered DGG systems plus CopernicusDEM(90), with a small level for conformance runs."
dgg_systems() = [
    ("HEALPix", DGG.HEALPixSystem(), 3),
    ("IGeo7", DGG.IGeo7System(), 3),
    ("H3", DGG.H3System(), 2),
    ("ISEA4R", DGG.ISEA4RSystem(), 3),
    ("S2", DGG.S2System(), 3),
    ("A5", DGG.A5System(), 2),
    ("CopernicusDEM(90)", DGG.CopernicusDEMSystem(90), 0),
]

end
