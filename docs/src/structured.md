```@meta
CurrentModule = GridInterface.Structured
```

# Structured grids

```@docs
Structured
```

## The face protocol

[`AbstractLatticeGrid`](@ref) supplies rectangular storage faces, indexing,
quadtrees and tiles. Its child [`AbstractStructuredGrid`](@ref) adds shared
lattice vertices, derived cell boundaries, analytic neighbours and vertex keys.
Both require [`facesize`](@ref) to return positive cell dimensions `(ni, nj)`.

A single-face shared-vertex provider implements three methods:

```julia
using GridInterface
import GridInterface.Structured as S
import GeometryOpsCore as GOC

struct SharedVertices{P<:AbstractMatrix} <: S.AbstractStructuredGrid
    points::P
end
GOC.best_manifold(::SharedVertices) = GOC.Planar()
S.facesize(g::SharedVertices, f) = (size(g.points, 1) - 1, size(g.points, 2) - 1)
S.facevertex(g::SharedVertices, f, a, b) = g.points[a, b]
```

The default [`nfaces`](@ref) is `1`. Default [`faceposition`](@ref) and
[`facecell`](@ref) use column-major positions (`i` fastest) within each face
and concatenate faces in face order. Counting `best_manifold`, the single-face
recipe needs three methods instead of six.

Several faces require a fourth method, `S.nfaces(g::MyGrid) = ...`, with
`facesize` and `facevertex` handling each face. The generic maps scan face sizes
and can be slower than provider-specific arithmetic. Explicit mapping methods
win dispatch; the concrete grids retain their optimized methods. Custom numbering
implements both maps and declares [`facecolumnmajor`](@ref) `false` when positions
do not ascend with `i` fastest. Seam and axis hooks remain optional.

A planar provider storing whole polygons uses the storage parent and three
methods, without shared-vertex geometry:

```julia
struct StoredPolygons{P<:AbstractMatrix} <: S.AbstractLatticeGrid
    polygons::P
end
GOC.best_manifold(::StoredPolygons) = GOC.Planar()
S.facesize(g::StoredPolygons, f) = size(g.polygons)
GridInterface.cell_geometry(g::StoredPolygons, k::Integer) = g.polygons[k]
```

It keeps the ordinary [`GeneralCells`](@ref GridInterface.GeneralCells) geometry
contract. Load GeometryOps to enable geometric extent, location and neighbour
fallbacks. Planar [`blockextent`](@ref) unions every cell extent in the block;
storage order does not imply adjacent polygon coordinates. Curved storage
providers must implement `blockextent(g, f, ir, jr)` with conservative bounds
covering every cell in the block.

These inherited defaults require the corresponding subtype. A foreign type with
another supertype needs an adapter; declaring [`isstructured`](@ref) alone does
not install the face, tree or tile methods.

```@docs
AbstractLatticeGrid
AbstractStructuredGrid
isstructured
nfaces
facesize
facevertex
faceposition
facecell
hasfaceaxes
faceaxes
facecoords
faceperiod
facelocalindex
facecolumnmajor
facecross
facefan
blockextent
```

## Seams

```@docs
FaceSeams
seamcross
```

## Trees and tiles

A single storage face has a quadtree; several faces add a root over their
quadtrees. Each node is an index block of a face. Its
[`position_span`](@ref GridInterface.position_span) is a
[`LatticeBlock`](@ref GridInterface.LatticeBlock) when positions are column-major,
so [`subtree`](@ref GridInterface.subtree) prunes it exactly against a tile.
Tiles are lazy `Subset`s whose positions stay in their type ([`TileIndices`](@ref)).
A tile's extent and its tree nodes use [`blockextent`](@ref); the cost of computing
bounds depends on the provider's implementation.

```@docs
TopDownQuadtreeCursor
FaceRootTree
LeafEntries
block_indices
TileIndices
TileChunks
```

## Concrete grids

```@docs
CellBasedGrid
RegularGrid
ExplicitPolygonGrid
```

## Internals

```@docs
Corners
circle_from_four_corners
```
