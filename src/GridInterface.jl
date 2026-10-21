"""
    GridInterface

One contract a grid implements once to be regridded, stencilled, analysed and
plotted by every consumer. The contract is a set of semi-orthogonal API
families (identity, tree, subsets and chunks, cell geometry and point location,
neighbours, dual grid, hierarchical systems); a grid opts into each.

**The required set is four methods:**

1. [`ncells`](@ref)`(grid)`;
2. [`best_manifold`](@ref)`(grid)` (`Planar()`, `Spherical(radius = r)`);
3. [`cellshape`](@ref)`(grid)`: `ConvexCells()` or `SimpleCells()` for
   single-ring grids (the default, `GeneralCells()`, refuses `cell_boundary`);
4. [`cell_boundary`](@ref)`(grid, i)` (single-ring grids) or
   [`cell_geometry`](@ref)`(grid, i)` (general grids).

With those four (and GeometryOps loaded) a grid gets centroids, areas,
extents, a tree, point location, neighbours, vertices, the dual, queries and
conformance. On a spherical (or ellipsoidal) manifold points must be
`UnitSphericalPoint`s and extents `SphericalCap`s (from
`UnitSphericalGeometryCore`, re-exported here); the core refuses other types
up front.

The contract is duck-typed: [`AbstractGrid`](@ref) is an optional convenience
supertype; types that already have a supertype elsewhere (DGG grids, CR
curvilinear grids) join by defining methods and `isgrid(::T) = true`.
GridInterface's own transparent wrappers ([`GridAndTree`](@ref),
[`ChunkedGrid`](@ref), [`GridAndVertices`](@ref)) share
[`AbstractWrappedGrid`](@ref) and forward every verb to their parent.

Positions (`1:ncells(grid)`) are the currency of every verb. Fallbacks that need
real geometry live in `GridInterfaceGeometryOpsExt` and activate when
GeometryOps is loaded; the core reaches them only through the hooks in
`fallbacks.jl`.
"""
module GridInterface

# --- The tree and unit-sphere imports ------------------------------------------
# SpatialTreeInterface owns the tree verbs (`split_weight`, `subtree_indices`
# included; GridInterface adds methods for its own trees and wrappers) and
# UnitSphericalGeometryCore the spherical point and cap types. Both are the
# sub-packages GeometryOps itself is built on, so nothing under src/ touches
# GeometryOps.
import SpatialTreeInterface as STI
import SpatialTreeInterface: split_weight, subtree_indices
import UnitSphericalGeometryCore: UnitSphericalPoint, SphericalCap, UnitSphereFromGeographic
# -------------------------------------------------------------------------------

import GeometryOpsCore as GOC
import GeometryOpsCore: Manifold, Planar, Spherical, best_manifold
import Extents
import GeoInterface as GI
# Small answers (neighbour rings, vertex fans, the corners of a cell with a
# bounded vertex count) are `SmallVector{N,T}`: capacity N, a runtime length,
# no heap storage. Re-exported so providers need no separate dependency.
using SmallCollections: SmallVector, push

include("fallbacks.jl")
include("lazy.jl")
include("identity.jl")
include("geometry.jl")
include("tree.jl")
include("subsets.jl")
include("restrict.jl")
include("neighbors.jl")
include("dual.jl")
include("hierarchy.jl")
include("cursor.jl")
include("multilevel.jl")
include("wrappers.jl")
include("structured/Structured.jl")

# --- identity
export AbstractGrid, isgrid
export AbstractCellIndex, iscellindex, ncells, cellindex, cellindextypes, reindex, localindex, globalindex, level
export Index, Local, Global, OUTSIDE
export SmallVector
# --- geometry and location
export CellShape, Polygonal, SingleRing, ConvexCells, SimpleCells, GeneralCells, cellshape
export EdgeKind, GeodesicEdges, CurvedEdges, edgekind
export cell_geometry, cell_boundary, cell_boundary!, cell_corners, cell_centroid, cell_area, cell_extent, best_manifold
export cellsize, cell_geometry_is_expensive, isglobal
export Tessellation, Tessellating, Overlapping, tessellation
export cellat, cellsat, cell_contains, location_is_analytic
export UnitSphericalPoint, SphericalCap
# --- tree
export treeify, subtree, split_weight, subtree_indices, position_span, span_range, span_mayhold, span_inside, LatticeBlock, leaf_capacity, grid_extent
export LeafHint, Balanced, ForSweeps, ForQueries
export query, query_candidates, query_region, QueryPredicate, Intersects, Within, Contains, CentroidIn
export RootIndex, root_index
# --- subsets and chunks
export subset, Subset, parentpositions, gridrelation, GridRelation, NoParent, ParentGrid, ParentSystem
export chunks, chunkat, grow, system, LevelChunks
# --- wrappers
export AbstractWrappedGrid, GridAndTree, ChunkedGrid, mockchunks, GridAndVertices
# --- neighbours
export Connectivity, Edge, Vertex, neighbors, ring, maxneighbors, maxring, winding, neighbors_are_analytic
export Winding, CounterClockwise, Clockwise, CustomOrder, Unordered
export halo, border, interior
# --- dual
export nvertices, cell_vertices, vertex_cells, vertex_point, vertex_isopen, dual, DualGrid, dual_ring_cells
export VertexKey, vertex_key, vertex_index, cell_vertex_keys, vertices_are_analytic, vertex_numbering_is_analytic, vertexat
export Conformity, Conforming, NonConforming, conformity
# --- hierarchy
export AbstractHierarchicalGridSystem, issystem, levels, rootcells, children, levelgrid
export has_sorted_subtrees, descendant_range, node_extent, has_congruent_refinement, ancestor
export HierarchyCursor, SelectionCursor, MaskedCursor, MultiLevelSet, hierarchy_chunks, expand_to_reference, aggregation_pairs

end
