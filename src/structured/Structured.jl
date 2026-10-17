"""
    GridInterface.Structured

The logically rectangular family: a grid that is one or more faces of an
(i, j) lattice of shared vertices (lon/lat and curvilinear rasters, cubed
spheres, HEALPix faces). A provider implements five hooks on its own type (a
subtype of [`AbstractStructuredGrid`](@ref)) and gets every GridInterface verb
derived: rings and corners, extents, a quadtree per face under a root over
faces, tile chunks, point location, aligned neighbours within and across
faces, analytic vertex keys.

    nfaces(grid); facesize(grid, f) -> (ni, nj); facevertex(grid, f, a, b) -> point
    faceposition(grid, f, i, j) -> position; facecell(grid, p) -> (f, i, j)

Faces with separable coordinates (a raster, lon/lat) add the axis hooks
([`hasfaceaxes`](@ref), [`faceaxes`](@ref), [`facecoords`](@ref),
[`faceperiod`](@ref), [`facelocalindex`](@ref)) and get O(1) point location.
Seams (a periodic axis, the edges of a cubed sphere) are declared by
arithmetic through [`facecross`](@ref); [`facefan`](@ref) covers the
vertices the walk across seams cannot reach (a lon/lat pole). A provider
with a closed-form block extent overrides [`blockextent`](@ref), which every
derived extent (tree nodes, leaf cells, tiles) reads.

The three concrete grids ([`CellBasedGrid`](@ref), [`RegularGrid`](@ref),
[`ExplicitPolygonGrid`](@ref)) and the quadtree
([`TopDownQuadtreeCursor`](@ref)) are ConservativeRegridding's, moved here.
"""
module Structured

import ..GridInterface as GrI
using ..GridInterface: AbstractGrid, ncells, cellshape, ConvexCells, GeneralCells, SingleRing,
    cell_boundary, cell_corners, cell_geometry, cell_extent, cell_contains, cellat, location_is_analytic,
    treeify, subtree, leaf_capacity, Balanced, ForSweeps, ForQueries, extent_is_expensive, grid_extent,
    chunks, chunkat, Subset, parentpositions,
    neighbors, Connectivity, Edge, Vertex, winding, CounterClockwise, Unordered, maxneighbors, neighbors_are_analytic,
    VertexKey, vertex_key, vertex_cells, vertices_are_analytic, isglobal,
    check_point, iscurved, UnitSphericalPoint, SphericalCap, UnitSphereFromGeographic
import ..GridInterface: STI, split_weight, subtree_indices
using SmallCollections: SmallVector, push
import GeometryOpsCore as GOC
using GeometryOpsCore: Manifold, Planar, Spherical, best_manifold
using UnitSphericalGeometryCore: GeographicFromUnitSphere, spherical_distance
import Extents
import GeoInterface as GI

include("faces.jl")
include("seams.jl")
include("grids.jl")
include("quadtree.jl")
include("tiles.jl")

export AbstractStructuredGrid, isstructured
export nfaces, facesize, facevertex, faceposition, facecell
export hasfaceaxes, faceaxes, facecoords, faceperiod, facelocalindex, facecolumnmajor
export facecross, facefan
export blockextent
export FaceSeams, seamcross, MaskedCursor
export CellBasedGrid, RegularGrid, ExplicitPolygonGrid

end
