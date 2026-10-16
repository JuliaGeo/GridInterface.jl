# CRAdapter, two directions, CR not edited:
#
# (a) CR's grids satisfy GridInterface: Trees.CellBasedGrid, Trees.RegularGrid,
#     Trees.ExplicitPolygonGrid (lon/lat and curvilinear quadtree grids) and, when
#     Healpix.jl is loaded, Healpix.HealpixMap through CR's Healpix extension.
#     Their storage partition is a lazy vector of raster tiles, each a typed
#     subset (`TileIndices`) whose extent, subtree and tree are O(1).
# (b) The bridge back: the unit CR's engine consumes is a grid paired with a
#     tree. `cr_tree(GridAndTree(grid, tree))` (or `cr_tree(grid)`, which pairs
#     the grid with `treeify(grid)`) is a `Trees`-compatible STI tree, so CR's
#     unmodified `Regridder` accepts ANY GridInterface grid:
#         Regridder(cr_tree(dst), cr_tree(src))            # any GI grid
#         Regridder(dst, src)                              # GI-native grid types (Subset, MultiLevelSet, GridAndTree, ...)
#         regridder(dst, src)                              # the same, refusing unsupported cell shapes up front
#     `split_weight` forwards through the bridge (GridInterface's for GI trees,
#     CR's own for CR/DGG trees), so CR's task split sees the real work estimate;
#     so does `subtree_indices`.
#
#   include("prototype/adapters/CRAdapter.jl"); using .CRAdapter
module CRAdapter

import GridInterface as GrI
using GridInterface
import ConservativeRegridding as CR
import ConservativeRegridding.Trees as Trees
import GeometryOps as GO
import GeometryOpsCore as GOC
using GeometryOps.UnitSpherical: UnitSphericalPoint
import GeoInterface as GI
using StaticArrays: SVector
import SpatialTreeInterface as STI

export CR, Trees, CRTree, cr_tree, regridder, lonlat_grid, cr_grid_names, TileIndices, TileChunks, TileCursor

# =============================================================================
# (b) the bridge: a Trees-compatible tree over a grid paired with a tree
# =============================================================================
"""
    CRTree(gt::GridAndTree)
    CRTree(gt, node)

STI tree that forwards traversal to `node` (the pair's tree, or a node of it,
whose leaves emit positions of the grid) and answers CR's cell-access hooks from
the grid: `getcell(t, i) = cell_geometry(grid, i)`, `cell_index_count(t) =
ncells(grid)`, `best_manifold(t) = best_manifold(grid)`,
`split_weight(t) = split_weight(node)`, `subtree_indices(t) = subtree_indices(node)`.
"""
struct CRTree{G,T}
    grid::G
    node::T
end
CRTree(gt::GridAndTree) = CRTree(gt, gt.tree)
cr_tree(gt::GridAndTree) = CRTree(gt)
cr_tree(g) = CRTree(GridAndTree(g))

STI.isspatialtree(::Type{<:CRTree}) = true
STI.node_extent_is_expensive(::Type{CRTree{G,T}}) where {G,T} = STI.node_extent_is_expensive(T)
STI.isleaf(t::CRTree) = STI.isleaf(t.node)
STI.nchild(t::CRTree) = STI.nchild(t.node)
STI.getchild(t::CRTree) = (CRTree(t.grid, c) for c in STI.getchild(t.node))
STI.getchild(t::CRTree, i::Int) = CRTree(t.grid, STI.getchild(t.node, i))
STI.node_extent(t::CRTree) = STI.node_extent(t.node)
STI.child_indices_extents(t::CRTree) = STI.child_indices_extents(t.node)

Trees.getcell(t::CRTree, i::Int) = GrI.cell_geometry(t.grid, i)
Trees.getcell(t::CRTree) = (GrI.cell_geometry(t.grid, i) for i in 1:GrI.ncells(t.grid))
# GridInterface trees answer GrI.split_weight, CR/DGG trees CR's Trees.split_weight
# (whose default is 1 on a GI tree, and GI's default is 1 on a foreign one).
_sw(n) = max(GrI.split_weight(n), Trees.split_weight(n))
Trees.split_weight(t::CRTree) = _sw(t.node)
Trees.ncells(t::CRTree) = _sw(t.node)
GrI.split_weight(t::CRTree) = _sw(t.node)
GrI.subtree_indices(t::CRTree) = GrI.subtree_indices(t.node)
Trees.cell_index_count(t::CRTree) = GrI.ncells(t.grid)
GOC.best_manifold(t::CRTree) = GOC.best_manifold(t.grid)

# GridInterface-native grid types (including GridAndTree) go straight into Regridder(dst, src).
Trees.treeify(::GOC.Manifold, g::GrI.AbstractGrid; kw...) = cr_tree(g)

"""
    regridder(dst, src; kw...) -> CR.Regridder

`CR.Regridder(cr_tree(dst), cr_tree(src); kw...)`, refusing up front what CR's
spherical engine cannot clip: it needs single-ring destination cells (it clips
source polygons against each destination ring, so convexity is destination-side).
`ConvexCells` is safe; `SimpleCells` destinations are accepted (DGG's densified
rings are convex in practice) and `GeneralCells` destinations are refused.
"""
function regridder(dst, src; kw...)
    m = GOC.best_manifold(dst)
    if m isa GOC.Spherical && !(GrI.cellshape(dst) isa GrI.SingleRing)
        throw(ArgumentError("CR's spherical engine clips against destination rings: " *
            "$(typeof(dst)) declares $(GrI.cellshape(dst)). Regrid onto a single-ring grid, " *
            "or use a planar view of the grid"))
    end
    return CR.Regridder(cr_tree(dst), cr_tree(src); kw...)
end

# =============================================================================
# (a) CR's curvilinear grids as GridInterface grids
# =============================================================================
const ACG = Trees.AbstractCurvilinearGrid
GrI.isgrid(::Type{<:ACG}) = true
_n1(g) = Trees.ncells(g, 1)
GrI.ncells(g::ACG) = Trees.ncells(g, 1) * Trees.ncells(g, 2)
@inline function _ij(g, k::Integer)
    1 <= k <= GrI.ncells(g) || throw(BoundsError(g, k))
    n = _n1(g)
    return (k - 1) % n + 1, (k - 1) ÷ n + 1
end
# Positions are CR's column-major linear index i + (j-1)·n1, which is what CR's
# quadtree cursor emits, so the cursor IS this grid's tree.
GrI.treeify(g::ACG) = Trees.TopDownQuadtreeCursor(g)
GrI.split_weight(q::Trees.TopDownQuadtreeCursor) = Trees.split_weight(q)
GrI.cell_extent(g::ACG, k::Integer) = ((i, j) = _ij(g, k); Trees.cell_range_extent(g, i:i, j:j))
GrI.extent_is_expensive(g::ACG) = Trees.extent_is_expensive(g)
GrI.cell_geometry_is_expensive(::ACG) = false

# geometry: corner points as stored (UnitSphericalPoint on the sphere)
GrI.cellshape(::Trees.CellBasedGrid) = GrI.ConvexCells()
GrI.cellshape(::Trees.RegularGrid) = GrI.ConvexCells()
GrI.cellshape(::Trees.ExplicitPolygonGrid) = GrI.GeneralCells()
function GrI.cell_boundary(g::Trees.CellBasedGrid, k::Integer)
    i, j = _ij(g, k); P = g.points
    return SVector(P[i, j], P[i+1, j], P[i+1, j+1], P[i, j+1])
end
_pt(::GOC.Planar, x, y) = (Float64(x), Float64(y))
_pt(::GOC.Spherical, x, y) = GO.UnitSphereFromGeographic()((x, y))
function GrI.cell_boundary(g::Trees.RegularGrid, k::Integer)
    i, j = _ij(g, k); m = GOC.manifold(g)
    return SVector(_pt(m, g.x[i], g.y[j]), _pt(m, g.x[i+1], g.y[j]), _pt(m, g.x[i+1], g.y[j+1]), _pt(m, g.x[i], g.y[j+1]))
end
GrI.cell_geometry(g::Trees.ExplicitPolygonGrid, k::Integer) = ((i, j) = _ij(g, k); Trees.getcell(g, i, j))
# a spherical lon/lat RegularGrid spanning 360° × 180° covers the sphere
GrI.isglobal(g::Trees.RegularGrid) = GOC.manifold(g) isa GOC.Spherical &&
    last(g.x) - first(g.x) >= 360 && first(g.y) <= -90 && last(g.y) >= 90

# planar RegularGrid: analytic location by binary search; ties to the lowest position
GrI.location_is_analytic(g::Trees.RegularGrid) = GOC.manifold(g) isa GOC.Planar
function GrI.cellat(g::Trees.RegularGrid, p)
    m = GOC.manifold(g)
    m isa GOC.Planar || return GrI.fallback_cellat(m, g, GrI.check_point(m, p))
    i = _bin(g.x, Float64(GI.x(p))); j = _bin(g.y, Float64(GI.y(p)))
    (i === nothing || j === nothing) && return nothing
    return i + (j - 1) * _n1(g)
end
function _bin(x, v)
    (first(x) <= v <= last(x)) || return nothing
    k = searchsortedfirst(x, v)
    return max(k - 1, 1)
end

# =============================================================================
# (a'') raster tiles: the storage partition as typed subsets
# =============================================================================
"""
    TileIndices(n1, ir, jr)

The column-major positions of the `ir × jr` tile of a grid with `n1` columns:
lazy and ascending, so `Subset(grid, TileIndices(...))` is a valid subset that
keeps its structure in its type. Its extent (`cell_range_extent`), its
`subtree` in the grid (CR's quadtree cursor over the tile, parent positions)
and its own tree (the same cursor with leaves remapped to local positions) are
all O(1).
"""
struct TileIndices <: AbstractVector{Int}
    n1::Int
    ir::UnitRange{Int}
    jr::UnitRange{Int}
end
Base.size(t::TileIndices) = (length(t.ir) * length(t.jr),)
Base.IndexStyle(::Type{TileIndices}) = IndexLinear()
@inline function Base.getindex(t::TileIndices, k::Int)
    @boundscheck checkbounds(t, k)
    a, b = divrem(k - 1, length(t.ir))
    return (first(t.ir) + b) + (first(t.jr) + a - 1) * t.n1
end
GrI._strictly_ascending(::TileIndices) = true
function Base.in(p::Integer, t::TileIndices)
    i = (p - 1) % t.n1 + 1; j = (p - 1) ÷ t.n1 + 1
    return p >= 1 && i in t.ir && j in t.jr
end

const TileSubset = Subset{<:ACG,TileIndices}

"""
    TileChunks(grid, (ti, tj))

Lazy vector of the `ti × tj` raster tiles of a CR curvilinear grid (edge tiles
smaller), column-major over tiles. `chunks(grid)` is `TileChunks(grid, (64, 64))`
(4096 cells, GR's default); `chunkat` is O(1).
"""
struct TileChunks{G} <: AbstractVector{Subset{G,TileIndices,Nothing}}
    grid::G
    ti::Int
    tj::Int
    nti::Int
    ntj::Int
end
function TileChunks(g::ACG, tile::NTuple{2,Int}=(64, 64))
    n1, n2 = Trees.ncells(g, 1), Trees.ncells(g, 2)
    return TileChunks(g, tile[1], tile[2], cld(n1, tile[1]), cld(n2, tile[2]))
end
Base.size(c::TileChunks) = (c.nti * c.ntj,)
Base.IndexStyle(::Type{<:TileChunks}) = IndexLinear()
function Base.getindex(c::TileChunks, k::Int)
    @boundscheck checkbounds(c, k)
    a, b = divrem(k - 1, c.nti)
    n1, n2 = Trees.ncells(c.grid, 1), Trees.ncells(c.grid, 2)
    ir = (b * c.ti + 1):min((b + 1) * c.ti, n1)
    jr = (a * c.tj + 1):min((a + 1) * c.tj, n2)
    return Subset(c.grid, TileIndices(n1, ir, jr))
end
function GrI.chunkat(c::TileChunks, p::Integer)
    i, j = _ij(c.grid, p)
    return (i - 1) ÷ c.ti + 1 + ((j - 1) ÷ c.tj) * c.nti
end
GrI.chunks(g::ACG) = TileChunks(g)
GrI.chunks(g::ACG, tile::NTuple{2,Int}) = TileChunks(g, tile)

GrI.grid_extent(s::TileSubset) = Trees.cell_range_extent(parent(s), parentpositions(s).ir, parentpositions(s).jr)
GrI.subtree(g::ACG, t::TileIndices) = Trees.TopDownQuadtreeCursor(g, (t.ir, t.jr))
GrI.subtree(g::ACG, s::TileSubset) = parent(s) === g ? GrI.subtree(g, parentpositions(s)) : GrI.subtree(g, collect(parentpositions(s)))

"Quadtree cursor over a tile whose leaves emit tile-local positions."
struct TileCursor{Q}
    q::Q
    ir0::Int
    jr0::Int
    ni::Int
    n1::Int
end
@inline function _loc(c::TileCursor, k::Int)
    i = (k - 1) % c.n1 + 1; j = (k - 1) ÷ c.n1 + 1
    return (i - c.ir0 + 1) + (j - c.jr0) * c.ni
end
_wrap(c::TileCursor, q) = TileCursor(q, c.ir0, c.jr0, c.ni, c.n1)
STI.isspatialtree(::Type{<:TileCursor}) = true
STI.node_extent_is_expensive(::Type{TileCursor{Q}}) where {Q} = STI.node_extent_is_expensive(Q)
STI.isleaf(c::TileCursor) = STI.isleaf(c.q)
STI.nchild(c::TileCursor) = STI.nchild(c.q)
STI.getchild(c::TileCursor) = (_wrap(c, q) for q in STI.getchild(c.q))
STI.getchild(c::TileCursor, i::Int) = _wrap(c, STI.getchild(c.q, i))
STI.node_extent(c::TileCursor) = STI.node_extent(c.q)
STI.child_indices_extents(c::TileCursor) = map(((k, e),) -> (_loc(c, k), e), STI.child_indices_extents(c.q))
GrI.split_weight(c::TileCursor) = GrI.split_weight(c.q)
GrI.subtree_indices(c::TileCursor) = map(k -> _loc(c, k), GrI.subtree_indices(c.q))
function GrI.treeify(s::TileSubset)
    t = parentpositions(s)
    return TileCursor(Trees.TopDownQuadtreeCursor(parent(s), (t.ir, t.jr)), first(t.ir), first(t.jr), length(t.ir), t.n1)
end

"""
    lonlat_grid(nlon, nlat; manifold = Spherical(), cellbased = true)

CR's standard global lon/lat test grid (corner matrix as in CR's test/trees/grids.jl).
"""
function lonlat_grid(nlon, nlat; manifold=GOC.Spherical(), cellbased=true)
    lons = range(-180, 180; length=nlon + 1); lats = range(-90, 90; length=nlat + 1)
    if cellbased
        pts = [GO.UnitSphereFromGeographic()((lo, la)) for lo in lons, la in lats]
        return Trees.CellBasedGrid(manifold, pts)
    end
    return Trees.RegularGrid(manifold, collect(lons), collect(lats))
end

cr_grid_names() = ["CellBasedGrid lon/lat (Spherical)", "RegularGrid (Planar)", "RegularGrid (Spherical)",
    "ExplicitPolygonGrid (Planar)", "CellBasedGrid warped curvilinear (Spherical)", "Healpix.HealpixMap (CR ext)"]

# =============================================================================
# (a') Healpix.HealpixMap through CR's Healpix extension (when Healpix is loaded)
# =============================================================================
"""
    adapt_healpix!(Healpix)

Load-time hook: call after `import Healpix` to make `Healpix.HealpixMap` a
GridInterface grid through CR's Healpix extension (its tree is CR's
HealpixRootNode; positions are the map's own pixel order + 1).
"""
function adapt_healpix!(Healpix::Module)
    @eval begin
        const HM = $(Healpix).HealpixMap
        GrI.isgrid(::Type{<:HM}) = true
        GrI.ncells(m::HM) = length(m)
        GrI.treeify(m::HM) = Trees.treeify(GO.Spherical(), m)
        const HExt = Base.get_extension(CR, :ConservativeRegriddingHealpixExt)
        GrI.split_weight(n::Union{HExt.HealpixRootNode,HExt.HealpixTreeNode}) = Trees.split_weight(n)
        GrI.cellshape(::HM) = GrI.SimpleCells()
        GrI.edgekind(::HM) = GrI.CurvedEdges()
        GrI.isglobal(::HM) = true
        GrI.cell_geometry(m::HM, i::Integer) = Trees.getcell(GrI.treeify(m), Int(i))
        function GrI.cell_boundary(m::HM, i::Integer)
            ring = collect(GI.getpoint(GI.getexterior(GrI.cell_geometry(m, i))))
            return ring[1:end-1]
        end
    end
end
export adapt_healpix!

end
