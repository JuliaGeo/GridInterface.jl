"""
    SpatialTreeInterfaceShim

Stand-in for the SpatialTreeInterface.jl package the PRD asks to split out of
GeometryOps, and for the light unit-sphere package the PRD lists as a
prerequisite split. It is the ONLY place the GridInterface stack touches
GeometryOps from a hard dependency. It re-exposes:

- `STI`: `GeometryOps.SpatialTreeInterface`;
- `UnitSpherical`: `GeometryOps.UnitSpherical` (`UnitSphericalPoint`,
  `SphericalCap`), the point and extent types every spherical grid must use;
- [`split_weight`](@ref): the traversal hint, owned here so it moves with the
  STI split and every wrapper (GridInterface's `GridAndTree`, consumer bridges)
  forwards one function.

Why a separate package rather than `import GeometryOps` inside GridInterface:
a Julia package extension cannot be triggered by a package that is also a hard
dependency, so the GeometryOps-backed fallbacks could not live in `ext/` if the
core depended on GeometryOps directly.
"""
module SpatialTreeInterfaceShim

import GeometryOps

const STI = GeometryOps.SpatialTreeInterface
const UnitSpherical = GeometryOps.UnitSpherical

"""
    split_weight(node) -> Int

Tree-node hint: a cheap estimate of the leaves under `node`, used to balance
dual traversals (which side to open first, how to split threaded work).
Wrong answers cost time, never correctness. Default `1`. Every tree wrapper
forwards it to the node it wraps.
"""
split_weight(node) = 1

export STI, UnitSpherical, split_weight

end
