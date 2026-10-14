"""
    SpatialTreeInterfaceShim

Stand-in for the SpatialTreeInterface.jl package the PRD asks to split out of
GeometryOps. It is the ONLY place the GridInterface stack touches GeometryOps
from a hard dependency: it re-exposes `GeometryOps.SpatialTreeInterface` as
`STI` and nothing else.

Why a separate package rather than `import GeometryOps` inside GridInterface:
a Julia package extension cannot be triggered by a package that is also a hard
dependency, so the GeometryOps-backed fallbacks could not live in `ext/` if the
core depended on GeometryOps directly. When SpatialTreeInterface.jl exists,
GridInterface replaces this dependency with it and nothing else changes.
"""
module SpatialTreeInterfaceShim

import GeometryOps

const STI = GeometryOps.SpatialTreeInterface

export STI

end
