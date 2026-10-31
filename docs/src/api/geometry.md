```@meta
CurrentModule = GridInterface
```

# Geometry

A grid declares the shape of its cells once ([`cellshape`](@ref)) and how boundary points relate to true edges ([`edgekind`](@ref)); the geometric verbs follow from the required boundary or geometry, with fallbacks in the GeometryOps extension.

```@docs
CellShape
Polygonal
SingleRing
ConvexCells
SimpleCells
GeneralCells
cellshape
EdgeKind
GeodesicEdges
CurvedEdges
edgekind
best_manifold
cell_geometry
cell_boundary
cell_boundary!
cell_corners
cell_centroid
cell_area
cell_extent
cellsize
cell_geometry_is_expensive
isglobal
Tessellation
Tessellating
Overlapping
tessellation
```

## Point location

```@docs
cellat
cellsat
cell_contains
location_is_analytic
```
