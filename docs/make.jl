using GridInterface
using Documenter

DocMeta.setdocmeta!(GridInterface, :DocTestSetup, :(using GridInterface); recursive=true)

makedocs(;
    modules=[GridInterface, GridInterface.Structured],
    authors="GridInterface contributors, Rafael Schouten, Anshul Singhvi, Gael Forget, Milan Klower, Simone Silvestri",
    sitename="GridInterface.jl",
    format=Documenter.HTML(;
        canonical="https://JuliaGeo.github.io/GridInterface.jl",
        edit_link="main",
        assets=String[],
        size_threshold=400 * 2^10,
    ),
    pages=[
        "Home" => "index.md",
        "Trees and restriction" => "trees.md",
        "API" => [
            "api.md",
            "Identity" => "api/identity.md",
            "Geometry" => "api/geometry.md",
            "Trees" => "api/trees.md",
            "Subsets, chunks and wrappers" => "api/subsets.md",
            "Neighbours" => "api/neighbors.md",
            "Vertices and the dual" => "api/dual.md",
            "Hierarchical systems" => "api/hierarchy.md",
        ],
        "Structured grids" => "structured.md",
    ],
)

deploydocs(;
    repo="github.com/JuliaGeo/GridInterface.jl",
    devbranch="main",
)
