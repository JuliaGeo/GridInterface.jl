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

## Start with the capabilities you need

Every recipe below supplies [`ncells`](@ref) and [`best_manifold`](@ref)
(`Planar()` or `Spherical(radius = r)`). With GeometryOps loaded:

| Capability | Additional provider methods | Generic operations |
|---|---|---|
| Bounds only: **3 methods total** | Conservative [`cell_extent`](@ref) | Spatial tree, restriction and candidate queries |
| Exact point location: **4 total** | `cell_extent` and [`cell_contains`](@ref) | `cellat` and `cellsat`, in addition to the bounds operations |
| General polygon geometry: **3 total** | [`cell_geometry`](@ref) | Measurements, bounds, containment, exact queries and applicable geometric neighbours |
| Single-ring geometry: **4 total** | [`cellshape`](@ref) and [`cell_boundary`](@ref) | Geometry fallbacks using the declared ring representation |

`cellshape` defaults to `GeneralCells()`. A single-ring provider declares
`ConvexCells()` or `SimpleCells()` according to its geometry. Vertex and dual
operations have additional shape and topology requirements; see
[Vertices and the dual](api/dual.md).

A provider can supply native adjacency with one [`neighbors`](@ref) override.
Discs, rings and halos have generic implementations; native algorithms and
capacity hints remain optional optimizations. See [Neighbours](api/neighbors.md).

On a spherical manifold points are `UnitSphericalPoint`s and extents
`SphericalCap`s; the core refuses other types up front.

These counts assume [`AbstractGrid`](@ref) registration. The core contract is
duck-typed: a foreign type adds `isgrid(::T) = true`. The
[Structured family](structured.md) has its own subtype-based extension contract.

## Julia versions

GridInterface requires Julia 1.11. Its zero-allocation guarantees (the allocation-free sweeps, trees and point location the conformance suite checks) hold on Julia 1.12 and later; on 1.11 the same code allocates a few bytes, and the conformance suite records those laws as known-broken.
