using SearchableGeometries
using Documenter

DocMeta.setdocmeta!(SearchableGeometries, :DocTestSetup, :(using SearchableGeometries); recursive=true)

makedocs(;
    modules=[SearchableGeometries],
    authors="Christina Taylor <cgtaylor@boisestate.edu>, Emmanuel Kwame Ayanful <emmanuelayanful@u.boisestate.edu>",
    sitename="SearchableGeometries.jl",
    format=Documenter.HTML(;
        canonical="https://cgt3.github.io/SearchableGeometries.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/cgt3/SearchableGeometries.jl",
    devbranch="main",
)
