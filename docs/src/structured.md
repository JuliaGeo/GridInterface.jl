```@meta
CurrentModule = GridInterface.Structured
```

# Structured grids

```@docs
Structured
```

## The face protocol

A provider subtypes [`AbstractStructuredGrid`](@ref) and implements the five face hooks; the seam and axis hooks are optional.

```@docs
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

The tree of a structured grid is a quadtree per face under a root over faces. A node is an index block of a face, so its [`position_span`](@ref GridInterface.position_span) is a [`LatticeBlock`](@ref GridInterface.LatticeBlock) when the face's positions are column-major: [`subtree`](@ref GridInterface.subtree) prunes it exactly against a tile. Tiles are the storage partition: lazy `Subset`s whose positions stay in their type ([`TileIndices`](@ref)), so a tile's extent, its subtree in the grid and its own tree are O(1).

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
