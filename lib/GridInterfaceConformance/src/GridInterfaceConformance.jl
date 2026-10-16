"""
    GridInterfaceConformance

Executable laws of GridInterface, one suite per API family (the PRD's law
table). Suites are gated by what a grid claims; a family the grid does not
claim is recorded as `:skip`, not failed. A generic fallback answering where
the grid claims a provider method (an analytic `cellat`, analytic neighbours)
is detected and failed; every verb answered by a fallback is listed.

    rep = check_grid(grid; name = "my grid")        # report only, never throws
    test_grid(grid; name = "my grid")               # same, as @testset/@test
    rep = check_system(sys; name = "HEALPix")       # hierarchy + performance + laziness laws
    rep = check_laziness(small, large)              # no verb's bytes scale between two sizes of a grid
    print_report(rep); family_status(rep)
"""
module GridInterfaceConformance

using Test
using Printf
import Random
import GridInterface as GrI
using GridInterface: ncells, cellindex, cellindextypes, localindex, globalindex, reindex,
    cell_geometry, cell_boundary, cell_centroid, cell_area, cell_extent, cellshape, best_manifold,
    treeify, subtree, cellat, neighbors, ring, maxneighbors, winding, tessellation, conformity,
    SingleRing, ConvexCells, GeneralCells, Tessellating, Overlapping, Conforming, NonConforming,
    Edge, Vertex, CounterClockwise, Clockwise
const STI = GrI.STI
import GeometryOps as GO
using GeometryOps.UnitSpherical: UnitSphericalPoint, SphericalCap, spherical_distance
import GeometryOpsCore as GOC
using GeometryOpsCore: Planar, Spherical
import Extents
import GeoInterface as GI

export check_grid, test_grid, check_system, test_system, check_multilevel, check_laziness,
    print_report, family_status, ConformanceReport, fallback_verbs

include("report.jl")
include("helpers.jl")
include("families.jl")
include("hierarchy.jl")

end
