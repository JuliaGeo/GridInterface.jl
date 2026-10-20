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
        "API" => "api.md",
        "Structured grids" => "structured.md",
    ],
)

deploydocs(;
    repo="github.com/JuliaGeo/GridInterface.jl",
    devbranch="main",
)
