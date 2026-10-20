```@meta
CurrentModule = GridInterface
```

# GridInterface

A shared contract for grids on the plane and the sphere: identity, geometry, neighbours, hierarchy, spatial trees and subsets, with conformance laws in `GridInterfaceConformance`. Providers include ConservativeRegridding.jl and DiscreteGlobalGrids.jl.

- [API](api.md): the core verbs and traits.
- [Structured grids](structured.md): the face protocol, quadtrees and tiles.
