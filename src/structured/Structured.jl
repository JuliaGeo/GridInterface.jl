"""
    GridInterface.Structured

[`AbstractLatticeGrid`](@ref) supplies rectangular storage faces, indexing,
quadtrees and tiles. Its child [`AbstractStructuredGrid`](@ref) derives geometry
from shared lattice vertices: implement `best_manifold`, [`facesize`](@ref) and
[`facevertex`](@ref) for a single face. Dimensions are positive cell counts.

[`nfaces`](@ref) defaults to `1`; [`faceposition`](@ref) and [`facecell`](@ref)
default to face-major positions with `i` fastest. Several faces add `nfaces`;
the generic maps scan face sizes, so providers may retain faster explicit maps.
Existing explicit methods, including the concrete grids' optimized maps, take
precedence. Custom numbering supplies both maps and declares
[`facecolumnmajor`](@ref) `false` when needed.

Storage-only planar providers implement `best_manifold`, `facesize` and
[`cell_geometry`](@ref GridInterface.cell_geometry). They keep ordinary geometry
traits and fallbacks; a block's extent unions every cell extent. Curved storage
providers must implement conservative [`blockextent`](@ref). These defaults
require the corresponding subtype; foreign types must adopt it or use an
adapter.

Shared-vertex faces with separable coordinates add the axis hooks
([`hasfaceaxes`](@ref), [`faceaxes`](@ref), [`facecoords`](@ref),
[`faceperiod`](@ref), [`facelocalindex`](@ref)) to derive point location from the
axes. Seams are declared through [`facecross`](@ref); [`facefan`](@ref) covers
vertices the seam walk cannot reach, such as a lon/lat pole. Providers with
closed-form block bounds override `blockextent`, which the default cell, tree
and tile extent implementations read.

The three concrete grids ([`CellBasedGrid`](@ref), [`RegularGrid`](@ref),
[`ExplicitPolygonGrid`](@ref)) and the quadtree
([`TopDownQuadtreeCursor`](@ref)) are ConservativeRegridding's, moved here.
"""
module Structured

import ..GridInterface as GrI
using ..GridInterface: AbstractGrid, ncells, cellshape, ConvexCells, GeneralCells, SingleRing,
    cell_boundary, cell_corners, cell_geometry, cell_extent, cell_contains, cellat, location_is_analytic,
    treeify, subtree, leaf_capacity, Balanced, ForSweeps, ForQueries, grid_extent,
    chunks, chunkat, Subset, parentpositions,
    neighbors, Connectivity, Edge, Vertex, winding, CounterClockwise, Unordered, maxneighbors, neighbors_are_analytic,
    VertexKey, vertex_key, vertex_cells, vertices_are_analytic, isglobal,
    check_point, iscurved, UnitSphericalPoint, SphericalCap, UnitSphereFromGeographic
import ..GridInterface: STI, split_weight, subtree_indices, position_span, MaskedCursor
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

export AbstractLatticeGrid, AbstractStructuredGrid, isstructured
export nfaces, facesize, facevertex, faceposition, facecell
export hasfaceaxes, faceaxes, facecoords, faceperiod, facelocalindex, facecolumnmajor
export facecross, facefan
export blockextent
# A provider whose storage tiles are irregular (a DiskArrays raster) builds its
# own `chunks` from these, and keeps O(1) chunk extents and subtrees.
export TileIndices
export FaceSeams, seamcross, MaskedCursor
export CellBasedGrid, RegularGrid, ExplicitPolygonGrid

end
