using GridInterface
using Documenter

DocMeta.setdocmeta!(GridInterface, :DocTestSetup, :(using GridInterface); recursive=true)

makedocs(;
    modules=[GridInterface],
    authors="GridInterface contributors, Rafael Schouten, Anshul Singhvi, Gael Forget, Milan Klower, Simone Silvestri",
    sitename="GridInterface.jl",
    format=Documenter.HTML(;
        canonical="https://JuliaGeo.github.io/GridInterface.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/JuliaGeo/GridInterface.jl",
    devbranch="main",
)
