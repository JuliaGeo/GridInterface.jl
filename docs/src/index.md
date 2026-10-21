```@meta
CurrentModule = GridInterface
```

# GridInterface

One contract a grid implements once to be regridded, stencilled, analysed and plotted by every consumer, on the plane and on the sphere. The contract is a set of API families; a grid opts into each, and `GridInterfaceConformance` (in `lib/`) checks a provider against the laws of the families it claims.

| family | what it answers | page |
|---|---|---|
| identity | `ncells`, positions `1:ncells(grid)`, typed ids | [Identity](api/identity.md) |
| geometry and location | boundaries, centroids, areas, extents, `cellat` | [Geometry](api/geometry.md) |
| tree | a spatial tree over the cells, queries | [Trees and restriction](trees.md) |
| subsets and chunks | a grid that remembers its parent; a declared partition | [Subsets, chunks and wrappers](api/subsets.md) |
| neighbours | adjacency, rings, halos | [Neighbours](api/neighbors.md) |
| dual | vertex identity and the dual grid | [Vertices and the dual](api/dual.md) |
| hierarchy | multi-level systems, level grids, multi-level sets | [Hierarchical systems](api/hierarchy.md) |
| structured | the face protocol: rasters, cubed spheres, quadtrees, tiles | [Structured grids](structured.md) |

## The required set

A grid joins with four methods:

1. [`ncells`](@ref)`(grid)`;
2. [`best_manifold`](@ref)`(grid)`: `Planar()` or `Spherical(radius = r)`;
3. [`cellshape`](@ref)`(grid)`: `ConvexCells()` or `SimpleCells()` for single-ring grids (the default, `GeneralCells()`, refuses `cell_boundary`);
4. [`cell_boundary`](@ref)`(grid, i)` for single-ring grids, or [`cell_geometry`](@ref)`(grid, i)` for general ones.

With those four, and GeometryOps loaded, a grid gets centroids, areas, extents, a tree, point location, neighbours, vertices, the dual, queries and conformance from the generic fallbacks. On a spherical manifold points are `UnitSphericalPoint`s and extents `SphericalCap`s; the core refuses other types up front.

The contract is duck-typed: [`AbstractGrid`](@ref) is an optional supertype, and a type with a supertype elsewhere joins by defining the verbs and `isgrid(::T) = true`.

## Julia versions

GridInterface requires Julia 1.11. Its zero-allocation guarantees (the allocation-free sweeps, trees and point location the conformance suite checks) hold on Julia 1.12 and later; on 1.11 the same code allocates a few bytes, and the conformance suite records those laws as known-broken.
