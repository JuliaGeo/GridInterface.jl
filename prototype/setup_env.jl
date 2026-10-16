# Builds prototype/env: develops GridInterface and its conformance package from
# this repository, GeometryOps with its SpatialTreeInterface and
# UnitSphericalGeometryCore sub-packages from the split worktree, and
# ConservativeRegridding, DiscreteGlobalGrids and GlobalRegridding from the
# sibling checkouts (read-only; nothing in them is modified). Run once:
#   julia prototype/setup_env.jl
using Pkg
const HERE = @__DIR__
const ROOT = dirname(HERE)
const GEO = dirname(ROOT)                      # /home/anshul/geo/gridinterface
const DGG = joinpath(GEO, "DiscreteGlobalGrids.jl")
const GO_SPLIT = "/home/anshul/.julia/dev/geo/GeometryOps-split"
Pkg.activate(joinpath(HERE, "env"))
Pkg.develop([
    PackageSpec(path=joinpath(GO_SPLIT, "SpatialTreeInterface")),
    PackageSpec(path=joinpath(GO_SPLIT, "UnitSphericalGeometryCore")),
    PackageSpec(path=GO_SPLIT),
    PackageSpec(path=ROOT),
    PackageSpec(path=joinpath(ROOT, "lib", "GridInterfaceConformance")),
    PackageSpec(path=joinpath(GEO, "ConservativeRegridding.jl")),
    PackageSpec(path=joinpath(DGG, "lib", "GlobalRegridding")),
    PackageSpec(path=DGG),
])
Pkg.add(["BenchmarkTools", "GeometryOpsCore", "Extents", "GeoInterface", "Test",
         "Healpix", "SparseArrays", "LinearAlgebra", "StaticArrays", "Random"])
Pkg.instantiate()
@time using GridInterface, DiscreteGlobalGrids, ConservativeRegridding
println("LOADED")
