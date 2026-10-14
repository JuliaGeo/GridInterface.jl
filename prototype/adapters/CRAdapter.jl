# CRAdapter, two directions, CR not edited:
#
# (a) CR's grids satisfy GridInterface: Trees.CellBasedGrid, Trees.RegularGrid,
#     Trees.ExplicitPolygonGrid (lon/lat and curvilinear quadtree grids) and, when
#     Healpix.jl is loaded, Healpix.HealpixMap through CR's Healpix extension.
#     (RingGrids, Oceananigans and ClimaCore grids need heavy packages that are
#     not in the depot; not adapted.)
# (b) The bridge back: `CRTree(grid, node)` is a `Trees`-compatible STI tree, so
#     CR's unmodified `Regridder` accepts ANY GridInterface grid:
#         Regridder(cr_tree(dst), cr_tree(src))            # any GI grid
#         Regridder(dst, src)                              # GI-native grid types (Subset, MultiLevelSet, ...)
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

const STI = GrI.STI
export CR, Trees, CRTree, cr_tree, lonlat_grid, cr_grid_names

# =============================================================================
# (b) the bridge: a Trees-compatible tree over any GridInterface grid
# =============================================================================
"""
    CRTree(grid, node)

STI tree that forwards traversal to `node` (a tree from `GridInterface.treeify`
or `subtree`, whose leaves emit positions of `grid`) and answers CR's cell-access
hooks from the grid: `getcell(t, i) = cell_geometry(grid, i)`,
`cell_index_count(t) = ncells(grid)`, `best_manifold(t) = best_manifold(grid)`.
The tree-side reference to the grid exists only in this bridge, because CR's
engine takes trees; GridInterface trees carry none.
"""
struct CRTree{G,T}
    grid::G
    node::T
end
cr_tree(g) = CRTree(g, GrI.treeify(g))

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
Trees.ncells(t::CRTree) = GrI.split_weight(t.node)
Trees.split_weight(t::CRTree) = GrI.split_weight(t.node)
Trees.cell_index_count(t::CRTree) = GrI.ncells(t.grid)
GOC.best_manifold(t::CRTree) = GOC.best_manifold(t.grid)

# GridInterface-native grid types go straight into Regridder(dst, src).
Trees.treeify(::GOC.Manifold, g::GrI.AbstractGrid; kw...) = cr_tree(g)

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

# planar RegularGrid: analytic location by binary search; ties to the lowest position
GrI.location_is_analytic(g::Trees.RegularGrid) = GOC.manifold(g) isa GOC.Planar
function GrI.cellat(g::Trees.RegularGrid, p)
    GOC.manifold(g) isa GOC.Planar || return GrI.fallback_cellat(GOC.manifold(g), g, p)
    i = _bin(g.x, Float64(GI.x(p))); j = _bin(g.y, Float64(GI.y(p)))
    (i === nothing || j === nothing) && return nothing
    return i + (j - 1) * _n1(g)
end
function _bin(x, v)
    (first(x) <= v <= last(x)) || return nothing
    k = searchsortedfirst(x, v)
    return max(k - 1, 1)
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
    adapt_healpix!()

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
        GrI.cellshape(::HM) = GrI.SimpleCells()
        GrI.cell_geometry(m::HM, i::Integer) = Trees.getcell(GrI.treeify(m), Int(i))
        function GrI.cell_boundary(m::HM, i::Integer)
            ring = collect(GI.getpoint(GI.getexterior(GrI.cell_geometry(m, i))))
            return ring[1:end-1]
        end
    end
end
export adapt_healpix!

end
