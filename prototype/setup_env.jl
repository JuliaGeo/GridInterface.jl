# Builds prototype/env: develops GridInterface, its conformance package and the
# STI shim from this repository, and ConservativeRegridding, DiscreteGlobalGrids
# and GlobalRegridding from the sibling checkouts (read-only; nothing in them
# is modified). Run once:
#   julia prototype/setup_env.jl
using Pkg
const HERE = @__DIR__
const ROOT = dirname(HERE)
const GEO = dirname(ROOT)                      # /home/anshul/geo/gridinterface
const DGG = joinpath(GEO, "DiscreteGlobalGrids.jl")
Pkg.activate(joinpath(HERE, "env"))
Pkg.develop([
    PackageSpec(path=joinpath(ROOT, "lib", "SpatialTreeInterfaceShim")),
    PackageSpec(path=ROOT),
    PackageSpec(path=joinpath(ROOT, "lib", "GridInterfaceConformance")),
    PackageSpec(path=joinpath(GEO, "ConservativeRegridding.jl")),
    PackageSpec(path=joinpath(DGG, "lib", "GlobalRegridding")),
    PackageSpec(path=DGG),
])
Pkg.add(["BenchmarkTools", "GeometryOps", "GeometryOpsCore", "Extents", "GeoInterface", "Test",
         "Healpix", "SparseArrays", "LinearAlgebra", "StaticArrays", "Random"])
Pkg.instantiate()
@time using GridInterface, DiscreteGlobalGrids, ConservativeRegridding
println("LOADED")
