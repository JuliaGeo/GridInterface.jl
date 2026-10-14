"""
    GridInterface

One contract a grid implements once to be regridded, stencilled, analysed and
plotted by every consumer. The contract is a set of semi-orthogonal API
families (identity, tree, subsets and chunks, cell geometry and point location,
neighbours, dual grid, hierarchical systems); a grid opts into each.

The contract is duck-typed: a grid is any object that answers
[`ncells`](@ref), [`best_manifold`](@ref) and either [`cell_boundary`](@ref)
(single-ring grids) or [`cell_geometry`](@ref). [`AbstractGrid`](@ref) is an
optional convenience supertype; types that already have a supertype elsewhere
(DGG grids, CR curvilinear grids) join by defining methods and
`isgrid(::Type{T}) = true`.

Positions (`1:ncells(grid)`) are the currency of every verb. Fallbacks that need
real geometry live in `GridInterfaceGeometryOpsExt` and activate when
GeometryOps is loaded; the core reaches them only through the hooks in
`fallbacks.jl`.
"""
module GridInterface

# --- The one tree import --------------------------------------------------------
# `STI` is GeometryOps.SpatialTreeInterface, re-exposed by a one-file shim package
# that stands in for SpatialTreeInterface.jl. Nothing else from GeometryOps is
# imported anywhere under src/.
using SpatialTreeInterfaceShim: STI
# -------------------------------------------------------------------------------

import GeometryOpsCore as GOC
import GeometryOpsCore: Manifold, Planar, Spherical, best_manifold
import Extents
import GeoInterface as GI

include("fallbacks.jl")
include("identity.jl")
include("geometry.jl")
include("tree.jl")
include("subsets.jl")
include("neighbors.jl")
include("dual.jl")
include("hierarchy.jl")
include("cursor.jl")
include("multilevel.jl")
include("gridandtree.jl")

# --- identity
export AbstractGrid, isgrid
export AbstractCellIndex, ncells, cellindex, cellindextypes, reindex, localindex, globalindex, level
export Index, Local, Global
# --- geometry and location
export CellShape, Polygonal, SingleRing, ConvexCells, SimpleCells, GeneralCells, cellshape
export cell_geometry, cell_boundary, cell_corners, cell_centroid, cell_area, cell_extent, best_manifold
export Tessellation, Tessellating, Overlapping, tessellation
export cellat, location_is_analytic
# --- tree
export treeify, subtree, GridAndTree, split_weight, extent_is_expensive, query, grid_extent
# --- subsets and chunks
export subset, Subset, system, chunks, parentpositions, gridrelation, NoParent, ParentGrid, ParentSystem
# --- neighbours
export Connectivity, Edge, Vertex, neighbors, ring, maxneighbors, winding, neighbors_are_analytic
export Winding, CounterClockwise, Clockwise, CustomOrder, Unordered
export halo, border, interior
# --- dual
export nvertices, cell_vertices, vertex_cells, vertex_point, dual, DualGrid
export Conformity, Conforming, NonConforming, conformity
# --- hierarchy
export AbstractHierarchicalGridSystem, issystem, levels, rootcells, children, levelgrid
export has_sorted_subtrees, descendant_range, node_extent, has_congruent_refinement, ancestor
export HierarchyCursor, SelectionCursor, MultiLevelSet, hierarchy_chunks

end
