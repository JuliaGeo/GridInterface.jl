# DGGAdapter: DiscreteGlobalGrids level grids, partial grids, multi-order sets
# and hierarchical systems satisfy GridInterface. Every method here is defined
# from outside on DGG's own types; DGG is not edited.
#
#   include("prototype/adapters/DGGAdapter.jl"); using .DGGAdapter
#
# What maps to what:
#   DGG.AbstractGrid (HierarchicalLevelGrid, PartialGrid, AuthalicGrid, ...)  → a GridInterface grid
#   DGG.AbstractHierarchicalGridSystem                                        → a GridInterface system
#   DGG.AbstractCellIndex                                                     → typed ids (declared with iscellindex)
#   DGG.MultiOrderCellSet                                                     → GridInterface.MultiLevelSet (O(1) wrap)
#
# Trees: the provider's own tree wins. `treeify` of a DGG grid is DGG's cursor
# (HierarchicalGridCursor, CopernicusDEM's MemoBlockCursor, ...). GridInterface's
# generic hierarchy cursor stays reachable as `GrI.hierarchy_tree(g)` (wrap with
# `GridAndTree`); on CopernicusDEM it uses a cached `RootIndex` over the 64,800 roots.
#
# Neighbours: DGG's analytic neighbours, disc form and ring bounds. HEALPix level
# grids are re-ordered here so that Edge neighbour k lies across edge
# (corner k, corner k+1) and declare `CounterClockwise`; the other systems keep
# DGG's counter-clockwise order, which starts elsewhere, so they declare
# `CustomOrder` (the alignment law is part of `CounterClockwise`).
#
# Vertices: HEALPix level grids name their vertices analytically
# (`vertices_are_analytic`): corner k of a pixel and the cells around it come
# from the nested lattice, with no table.
module DGGAdapter

import GridInterface as GrI
using GridInterface
import DiscreteGlobalGrids as DGG
const FB = DGG.Fallbacks
const HP = DGG.HEALPix
import GeometryOpsCore as GOC
import ConservativeRegridding.Trees as Trees

export dgg_multilevel, dgg_tree, dgg_systems, DGG

const G = DGG.AbstractGrid
const SYS = DGG.AbstractHierarchicalGridSystem
const ID = DGG.AbstractCellIndex
const HLG = DGG.HierarchicalLevelGrid
const CDEM = DGG.CopernicusDEM.CopernicusDEMSystem
const HPG = HLG{HP.HEALPixSystem}

# --- recognition --------------------------------------------------------------
GrI.isgrid(::G) = true
GrI.issystem(::SYS) = true
GrI.iscellindex(::ID) = true      # DGG ids cannot subtype GrI.AbstractCellIndex
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
GrI.isglobal(::HLG) = true

# --- geometry -----------------------------------------------------------------
# Cell-shape class and edge kind, declared once per system (they hold for every
# level grid). ConvexCells/GeodesicEdges where the boundary is the great-circle
# polygon of the cell; systems whose true edges are chart curves densify them
# and declare SimpleCells/CurvedEdges.
dgg_cellshape(::Any) = GrI.SimpleCells()
dgg_cellshape(::DGG.S2System) = GrI.ConvexCells()
dgg_cellshape(::DGG.H3System) = GrI.ConvexCells()
dgg_cellshape(::DGG.IGeo7System) = GrI.ConvexCells()
# CopernicusDEM cells are lon/lat boxes: parallels are not great circles.
dgg_cellshape(s::DGG.AuthalicSystem) = dgg_cellshape(s.system)
dgg_edgekind(s) = dgg_cellshape(s) isa GrI.ConvexCells ? GrI.GeodesicEdges() : GrI.CurvedEdges()
dgg_edgekind(::DGG.H3System) = GrI.GeodesicEdges()       # bent at class III levels, but each piece is geodesic
dgg_edgekind(s::DGG.AuthalicSystem) = dgg_edgekind(s.system)
GrI.cellshape(g::G) = dgg_cellshape(DGG.system(g))
# H3 is convex only at even (class II) levels: class III cells cross icosahedron
# edges and bend there. So the shape is per level here, not per system.
GrI.cellshape(g::HLG{DGG.H3System}) = iseven(g.level) ? GrI.ConvexCells() : GrI.SimpleCells()
GrI.edgekind(g::G) = dgg_edgekind(DGG.system(g))
GrI.extent_is_expensive(g::G) = !(FB.cell_cap_is_cheap(g) isa Val{true})
# DGG builds a polygon per call (densified chart edges, H3/IGeo7 projections):
# worth memoising per chunk. CopernicusDEM cells are four-corner boxes.
GrI.cell_geometry_is_expensive(g::G) = !(DGG.system(g) isa CDEM)
GrI.cellsize(g::G) = DGG.cellsize(g; radius=Float64(GOC.best_manifold(g).radius))

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
GrI.cell_corners(g::G, c::ID) = DGG.cell_corners(g, c)
GrI.cell_centroid(g::G, c::ID) = DGG.cell_centroid(g, c)
GrI.cell_area(g::G, c::ID) = DGG.cell_area(g, c)
GrI.cell_extent(g::G, c::ID) = FB.cell_cap(g, c)

# --- location -----------------------------------------------------------------
function GrI.location_is_analytic(g::G)
    s = DGG.system(g)
    return s !== nothing && DGG.has_direct_location(s)
end
function GrI.cellat(g::G, p)
    m = GOC.best_manifold(g)
    GrI.location_is_analytic(g) || return GrI.fallback_cellat(m, g, GrI.check_point(m, p))
    c = DGG.cellat(g, GrI.check_point(m, p))
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
GrI.neighbors(g::G, i::Integer, k::Integer; connectivity::GrI.Connectivity=GrI.Vertex()) =
    DGG.neighbors(g, Int(i), k; connectivity=_dconn(connectivity))
GrI.ring(g::G, i::Integer, k::Integer; connectivity::GrI.Connectivity=GrI.Vertex()) =
    DGG.ring(g, Int(i), k; connectivity=_dconn(connectivity))
GrI.maxneighbors(g::G, c::GrI.Connectivity) = DGG.maxneighbors(g, _dconn(c))
GrI.maxring(g::G, k::Integer, c::GrI.Connectivity=GrI.Vertex()) =
    k == 0 ? 1 : k == 1 ? DGG.maxneighbors(g, _dconn(c)) : DGG.maxring(g, k, _dconn(c))
# DGG's counter-clockwise rings start on a system-specific spoke, not across
# edge (corner 1, corner 2): only the re-ordered HEALPix grids keep CounterClockwise.
GrI.winding(::HPG, ::GrI.Connectivity) = GrI.CounterClockwise()
function GrI.winding(g::G, c::GrI.Connectivity)
    w = _gwind(DGG.winding(g, _dconn(c)))
    return w isa GrI.CounterClockwise && !(g isa HPG) ? GrI.CustomOrder() : w
end

"A fixed-capacity vector of Int on the stack (neighbour rings, vertex fans)."
struct SmallInts{N} <: AbstractVector{Int}
    v::NTuple{N,Int}
    n::Int
end
Base.size(s::SmallInts) = (s.n,)
Base.IndexStyle(::Type{<:SmallInts}) = IndexLinear()
Base.@propagate_inbounds Base.getindex(s::SmallInts, i::Int) = (@boundscheck checkbounds(s, i); s.v[i])

# HEALPix: nested_neighbors gives compass slots SW, W, NW, N, NE, E, SE, S (1..8),
# -1 where a degree-3 vertex has no diagonal neighbour; level 0 repeats.
# Corners are N, W, S, E; edge k = (corner k, corner k+1).
const _HP_EDGE = (3, 1, 7, 5)                   # NW, SW, SE, NE: across edges NW, WS, SE, EN
const _HP_VERT = (3, 2, 1, 8, 7, 6, 5, 4)       # NW, W, SW, S, SE, E, NE, N
@inline function _hp_pick(nb, i::Int, slots::NTuple{K,Int}) where {K}
    out = ntuple(_ -> 0, Val(K)); n = 0
    for m in slots
        q = nb[m]
        q < 0 && continue
        p = Int(q) + 1
        p == i && continue
        dup = false
        for j in 1:n
            out[j] == p && (dup = true)
        end
        dup && continue
        n += 1
        out = Base.setindex(out, p, n)
    end
    return SmallInts(out, n)
end
@inline function _hp_nb(g::HPG, i::Integer)
    1 <= i <= DGG.ncells(g) || throw(BoundsError(g, i))
    return HP.nested_neighbors(Int(i) - 1, g.level)
end
function GrI.neighbors(g::HPG, i::Integer; connectivity::GrI.Connectivity=GrI.Vertex())
    nb = _hp_nb(g, i)
    return connectivity isa GrI.Edge ? SmallInts(_pad8(_hp_pick(nb, Int(i), _HP_EDGE))) : _hp_pick(nb, Int(i), _HP_VERT)
end
_pad8(s::SmallInts{4}) = (ntuple(k -> k <= 4 ? s.v[k] : 0, Val(8)), s.n)
SmallInts(t::Tuple{NTuple{8,Int},Int}) = SmallInts{8}(t[1], t[2])
function GrI.neighbors(g::HPG, i::Integer, k::Integer; connectivity::GrI.Connectivity=GrI.Vertex())
    r1 = GrI.neighbors(g, i; connectivity)
    k <= 1 && return k == 1 ? collect(r1) : Int[]
    d = DGG.neighbors(g, Int(i), k; connectivity=_dconn(connectivity))
    return vcat(collect(r1), d[length(r1)+1:end])
end
GrI.ring(g::HPG, i::Integer, k::Integer; connectivity::GrI.Connectivity=GrI.Vertex()) =
    k == 1 ? collect(GrI.neighbors(g, i; connectivity)) : DGG.ring(g, Int(i), k; connectivity=_dconn(connectivity))

# --- region verbs: DGG's O(border) engines on rooted subtrees ---------------------
# The primitive `halo(parent, positions, k)` on a level grid: when `positions` is
# the descendant range of one cell (a rooted chunk), DGG's subtree engine answers.
# (Not CopernicusDEM: DGG's subtree halo engine does not finish on a 1° tile of
# 1.06M cells; the generic merge takes seconds there.)
function _rooted(lg::HLG, r::AbstractUnitRange{<:Integer})
    sys = DGG.system(lg)
    (isempty(r) || !DGG.has_sorted_subtrees(sys) || sys isa CDEM) && return nothing
    c = DGG.cellindex(lg, Int(first(r)))
    for l in DGG.levels(sys)
        l > lg.level && break
        a = DGG.ancestor(sys, c, l)
        DGG.descendant_range(sys, a, lg.level) == r && return DGG.subtree(sys, a, lg.level)
    end
    return nothing
end
const _GENERIC_HALO = Tuple{Any,AbstractVector{<:Integer},Integer}
const _GENERIC_BORDER = Tuple{Any,AbstractVector{<:Integer}}
function GrI.halo(lg::HLG, positions::AbstractUnitRange{<:Integer}, k::Integer=1; connectivity::GrI.Connectivity=GrI.Vertex())
    pg = k == 1 && connectivity isa GrI.Vertex ? _rooted(lg, positions) : nothing
    pg === nothing && return invoke(GrI.halo, _GENERIC_HALO, lg, positions, k; connectivity)
    return sort!(collect(Int, DGG.halo(pg)))
end
function GrI.border(lg::HLG, positions::AbstractUnitRange{<:Integer}; connectivity::GrI.Connectivity=GrI.Vertex())
    pg = connectivity isa GrI.Vertex ? _rooted(lg, positions) : nothing
    pg === nothing && return invoke(GrI.border, _GENERIC_BORDER, lg, positions; connectivity)
    return sort!(collect(Int, DGG.border(pg)))
end
GrI.halo(pg::DGG.PartialGrid, k::Integer=1; connectivity::GrI.Connectivity=GrI.Vertex()) =
    k == 1 && connectivity isa GrI.Vertex ? sort!(collect(Int, DGG.halo(pg))) :
    GrI.halo(parent(pg), GrI.parentpositions(pg), k; connectivity)
GrI.border(pg::DGG.PartialGrid; connectivity::GrI.Connectivity=GrI.Vertex()) =
    connectivity isa GrI.Vertex ? sort!(collect(Int, DGG.border(pg))) :
    GrI.border(parent(pg), GrI.parentpositions(pg); connectivity)

# --- trees and chunks ---------------------------------------------------------------
# The provider's own tree wins over the generic hierarchy cursor. DGG's cursors
# have one shape, so the capacity every intent maps to is ignored (O(1) as before).
GrI.treeify(g::G, ::Integer) = DGG.treeify(g)
"DGG's own tree for a grid (also what `treeify` returns): `GridAndTree(g, dgg_tree(g))`."
dgg_tree(g) = DGG.treeify(g)
# DGG's cursors size themselves through CR's `Trees.ncells`; forward it.
const DGGTree = Union{DGG.Engine.HierarchicalGridCursor,DGG.Engine.IndexTreeNode,
    DGG.CopernicusDEM.BlockCursor,DGG.CopernicusDEM.MemoBlockCursor}
GrI.split_weight(t::DGGTree) = Trees.split_weight(t)

# CopernicusDEM has 64,800 roots: the generic cursor gets a packed root index, built once.
const _ROOT_INDEX = IdDict{Any,Any}()
const _ROOT_LOCK = ReentrantLock()
GrI.root_index(s::CDEM) = lock(() -> get!(() -> GrI.RootIndex(s), _ROOT_INDEX, s), _ROOT_LOCK)::GrI.RootIndex

"""
    chunks(g)

The storage partition of a DGG level or partial grid on a sorted system: the
hierarchy chunks at the level whose chunks hold closest to 4096 cells (GR's
default). CopernicusDEM: one chunk per 1° tile (level 0). `chunks(g, l)` picks
the level. A5 (unsorted) declares none.
"""
function GrI.chunks(g::G)
    s = DGG.system(g)
    (s === nothing || !DGG.has_sorted_subtrees(s)) && return nothing
    return GrI.hierarchy_chunks(g, dgg_chunk_level(s, DGG.level(g)))
end
function dgg_chunk_level(s, l)
    n = Float64(DGG.ncells(DGG.levelgrid(s, l)))
    best, bestd = l, Inf
    for lc in DGG.levels(s)
        lc > l && break
        d = abs(log(n / DGG.ncells(DGG.levelgrid(s, lc))) - log(4096))
        d < bestd && ((best, bestd) = (lc, d))
    end
    return best
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

# --- analytic vertices: HEALPix -------------------------------------------------------
# Around corner k (N, W, S, E) of pixel p, counter-clockwise from p: the compass
# slots of the three other cells (one is missing at a degree-3 vertex).
const _HP_FAN = ((5, 4, 3), (3, 2, 1), (1, 8, 7), (7, 6, 5))
GrI.vertices_are_analytic(::HPG) = true

@inline function _hp_corner(g::HPG, i::Int, k::Int)
    nside = 1 << g.level
    ix, iy, f = HP.nested_to_xyf(i - 1, nside)
    return HP.pixel_corners(ix, iy, f, nside)[k]
end
@inline function _hp_fan(g::HPG, i::Int, k::Int)
    1 <= k <= 4 || throw(BoundsError(1:4, k))
    nb = _hp_nb(g, i)
    s = _hp_pick(nb, i, _HP_FAN[k])
    return SmallInts((i, s.v[1], s.v[2], s.v[3]), s.n + 1)
end
@inline _d2(a, b) = (a[1] - b[1])^2 + (a[2] - b[2])^2 + (a[3] - b[3])^2
# lattice offset of corner k (N, W, S, E) from the pixel's (ix, iy)
const _HP_COFF = ((1, 1), (0, 1), (0, 0), (1, 0))
function GrI.vertex_key(g::HPG, i::Integer, k::Integer)
    fan = _hp_fan(g, Int(i), Int(k))
    o = minimum(fan)
    o == i && return GrI.VertexKey(Int(i), Int(k))
    nside = 1 << g.level
    ix, iy, f = HP.nested_to_xyf(Int(i) - 1, nside)
    ox, oy, fo = HP.nested_to_xyf(o - 1, nside)
    if fo == f          # same face: the corner index follows from the lattice
        vx, vy = ix + _HP_COFF[k][1], iy + _HP_COFF[k][2]
        for j in 1:4
            (ox + _HP_COFF[j][1], oy + _HP_COFF[j][2]) == (vx, vy) && return GrI.VertexKey(o, j)
        end
    end
    p = HP.pixel_corners(ix, iy, f, nside)[k]       # across a face edge: nearest corner
    cs = HP.pixel_corners(ox, oy, fo, nside)
    j = 1; d = _d2(cs[1], p)
    for jj in 2:4
        dj = _d2(cs[jj], p)
        dj < d && ((j, d) = (jj, dj))
    end
    return GrI.VertexKey(o, j)
end
GrI.vertex_cells(g::HPG, key::GrI.VertexKey) = _hp_fan(g, key.cell, key.corner)
GrI.vertex_point(g::HPG, key::GrI.VertexKey) = _hp_corner(g, key.cell, key.corner)

# Dense vertex numbers, O(1). Lattice point (f, x, y), x, y ∈ 0:N, N = nside. Every
# face owns x ∈ 0:N-1, y ∈ 1:N (its SW and NW edges and W corner, not its NE or SE
# edge): each shared edge is NE/SE on one side and SW/NW on the other (face tables
# in DGG's neighbors.jl), so each non-pole vertex has one owner. The poles (N corner
# of faces 0:3, S corner of 8:11) are left over and numbered last.
#   v = f N² + (y - 1) N + x + 1,   north pole 12N² + 1,   south pole 12N² + 2
GrI.vertex_numbering_is_analytic(::HPG) = true

# Owner of lattice point (f, x, y), or face -1 / -2 for the north / south pole.
# Crossing NE (x = N) or SE (y = 0): north f ↔ (f+1)%4 via (N,y) ↔ (y,N) and
# ↔ 4+(f+1)%4 via (x,0) ↔ (x,N); equatorial 4+j ↔ north j via (N,y) ↔ (0,y) and
# ↔ 8+j via (x,0) ↔ (x,N); south 8+j ↔ 4+(j+1)%4 via (N,y) ↔ (0,y) and
# ↔ 8+(j+1)%4 via (x,0) ↔ (0,x). An E corner (N, 0) crosses twice.
@inline function _hp_owner(f::Int, x::Int, y::Int, N::Int)
    r, j1 = f >> 2, (f + 1) & 3
    if x == N
        if r == 0
            y == N && return (-1, 0, 0)
            return (j1, y, N)
        end
        y > 0 && return (r == 1 ? f & 3 : 4 + j1, 0, y)
        return (r == 1 ? 4 + j1 : 8 + j1, 0, N)
    elseif y == 0
        r == 0 && return (4 + j1, x, N)
        r == 1 && return (f + 4, x, N)
        x == 0 && return (-2, 0, 0)
        return (8 + j1, 0, x)
    end
    return (f, x, y)
end
function GrI.nvertices(g::HPG)     # Euler: V = F + 2 for a quad tessellation of the sphere
    N = 1 << g.level
    return 12 * N * N + 2
end
function GrI.vertex_key(g::HPG, v::Integer)
    lv = g.level; N = 1 << lv; nf = N * N
    1 <= v <= 12nf + 2 || throw(BoundsError(1:12nf+2, v))
    v > 12nf && return v == 12nf + 1 ? GrI.vertex_key(g, nf, 1) : GrI.vertex_key(g, 8nf + 1, 3)
    w = Int(v) - 1
    f, r = w >> 2lv, w & (nf - 1)
    x, y = r & (N - 1), r >> lv + 1
    return GrI.vertex_key(g, HP.xyf_to_nested(x, y - 1, f, N) + 1, 2)   # owner point is the W corner of (x, y-1)
end
function GrI.vertex_index(g::HPG, key::GrI.VertexKey)
    lv = g.level; N = 1 << lv
    k = key.corner
    1 <= k <= 4 || throw(BoundsError(1:4, k))
    ix, iy, f = HP.nested_to_xyf(key.cell - 1, N)
    fo, x, y = _hp_owner(f, ix + _HP_COFF[k][1], iy + _HP_COFF[k][2], N)
    fo < 0 && return 12 * N * N - fo
    return (fo << 2lv) + ((y - 1) << lv) + x + 1
end

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

"The registered DGG systems plus CopernicusDEM(90), with a small level for conformance runs (H3 at an odd level too)."
dgg_systems() = [
    ("HEALPix", DGG.HEALPixSystem(), 3),
    ("IGeo7", DGG.IGeo7System(), 3),
    ("H3", DGG.H3System(), 2),
    ("H3 (odd level)", DGG.H3System(), 1),
    ("ISEA4R", DGG.ISEA4RSystem(), 3),
    ("S2", DGG.S2System(), 3),
    ("A5", DGG.A5System(), 2),
    ("CopernicusDEM(90)", DGG.CopernicusDEMSystem(90), 0),
]

end
