# GridInterface

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaGeo.github.io/GridInterface.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaGeo.github.io/GridInterface.jl/dev/)
[![Build Status](https://github.com/JuliaGeo/GridInterface.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JuliaGeo/GridInterface.jl/actions/workflows/CI.yml?query=branch%3Amain)

A shared contract for grids on the plane and the sphere: identity, geometry, neighbours, hierarchy, spatial trees and subsets. `lib/GridInterfaceConformance` checks a provider against the contract's laws.

GridInterface builds on the unregistered `SpatialTreeInterface` and `UnitSphericalGeometryCore` from GeometryOps.jl's `as/split-unitsphericalgeometrycore` branch, through `[sources]`, so it needs Julia 1.11 or later.
