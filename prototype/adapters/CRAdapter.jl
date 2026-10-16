# CRAdapter, two directions, CR not edited:
#
# (a) CR's grids become GridInterface grids by construction: a CR
#     `Trees.CellBasedGrid`, `Trees.RegularGrid` or `Trees.ExplicitPolygonGrid`
#     is converted to its `GridInterface.Structured` twin (`structured(g)`),
#     which shares its corners, axes or polygons and answers every verb from
#     the face hooks (quadtree, tiles, location, neighbours, vertices). When
#     Healpix.jl is loaded, Healpix.HealpixMap joins through CR's Healpix
#     extension.
# (b) The bridge back: the unit CR's engine consumes is a grid paired with a
#     tree. `cr_tree(GridAndTree(grid, tree))` (or `cr_tree(grid)`, which pairs
#     the grid with `treeify(grid; leaf_capacity = ForSweeps())`, the sweep
#     intent) is a `Trees`-compatible STI tree, so CR's
#     unmodified `Regridder` accepts ANY GridInterface grid:
#         Regridder(cr_tree(dst), cr_tree(src))            # any GI grid
#         Regridder(dst, src)                              # GI-native grid types (Subset, MultiLevelSet, GridAndTree, Structured grids, ...)
#         regridder(dst, src)                              # the same, refusing unsupported cell shapes up front
#     `split_weight` forwards through the bridge (GridInterface's for GI trees,
#     CR's own for CR/DGG trees), so CR's task split sees the real work estimate;
#     so does `subtree_indices`.
#
#   include("prototype/adapters/CRAdapter.jl"); using .CRAdapter
module CRAdapter

import GridInterface as GrI
using GridInterface
using GridInterface.Structured
import ConservativeRegridding as CR
import ConservativeRegridding.Trees as Trees
import GeometryOps as GO
import GeometryOpsCore as GOC
import GeoInterface as GI
import SpatialTreeInterface as STI

export CR, Trees, CRTree, cr_tree, regridder, structured, lonlat_grid, cr_lonlat_grid, cr_grid_names

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
cr_tree(g) = CRTree(GridAndTree(g; leaf_capacity=GrI.ForSweeps()))   # CR runs a dual-tree sweep

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

# GridInterface-native grid types (including the wrappers) go straight into Regridder(dst, src).
Trees.treeify(::GOC.Manifold, g::Union{GrI.AbstractGrid,GrI.AbstractWrappedGrid}; kw...) = cr_tree(g)

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
# (a) CR's curvilinear grids: their Structured twins, by construction
# =============================================================================
"""
    structured(g::Trees.AbstractCurvilinearGrid; periodic = nothing, collapsed = nothing)

The `GridInterface.Structured` twin of a CR grid: the same manifold and the
same corner matrix, axis vectors or polygon matrix (shared, not copied).
Positions are CR's column-major linear index `i + (j - 1) n1`. CR's grids
carry no topology, so the seams (`periodic`, `collapsed`: see
`Structured.CellBasedGrid`) are detected from the vertices unless given.
"""
structured(g::Trees.CellBasedGrid; kw...) = Structured.CellBasedGrid(g.manifold, g.points; kw...)
structured(g::Trees.RegularGrid; kw...) = Structured.RegularGrid(g.manifold, g.x, g.y; kw...)
structured(g::Trees.ExplicitPolygonGrid) = Structured.ExplicitPolygonGrid(g.manifold, g.polygons)
structured(g::GrI.AbstractGrid) = g

"""
    cr_lonlat_grid(nlon, nlat; manifold = Spherical(), cellbased = true)

CR's standard global lon/lat test grid (corner matrix as in CR's test/trees/grids.jl), as CR's own type.
"""
function cr_lonlat_grid(nlon, nlat; manifold=GOC.Spherical(), cellbased=true)
    lons = range(-180, 180; length=nlon + 1); lats = range(-90, 90; length=nlat + 1)
    if cellbased
        pts = [GO.UnitSphereFromGeographic()((lo, la)) for lo in lons, la in lats]
        return Trees.CellBasedGrid(manifold, pts)
    end
    return Trees.RegularGrid(manifold, collect(lons), collect(lats))
end

"""
    lonlat_grid(nlon, nlat; manifold = Spherical(), cellbased = true)

The same grid as a GridInterface grid: the `Structured` twin of [`cr_lonlat_grid`](@ref),
with its seams declared rather than detected: longitude wraps at the dateline
and both pole rows collapse (on a plane the grid is a regional raster).
"""
function lonlat_grid(nlon, nlat; manifold=GOC.Spherical(), cellbased=true)
    g = cr_lonlat_grid(nlon, nlat; manifold, cellbased)
    GrI.iscurved(manifold) || return structured(g; periodic=(false, false), collapsed=(false, false))
    return structured(g; periodic=(true, false), collapsed=(true, true))
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
        GrI.isgrid(::HM) = true
        GrI.ncells(m::HM) = length(m)
        GrI.treeify(m::HM, ::Integer) = Trees.treeify(GO.Spherical(), m)
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
