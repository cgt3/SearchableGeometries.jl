using SearchableGeometries
using Documenter

DocMeta.setdocmeta!(
    SearchableGeometries,
    :DocTestSetup,
    :(using SearchableGeometries);
    recursive=true
)

makedocs(
    modules=[SearchableGeometries],
    authors="Christina Taylor <cgtaylor@boisestate.edu>, Emmanuel Kwame Ayanful <emmanuelayanful@u.boisestate.edu>",
    sitename="SearchableGeometries.jl",
    format=Documenter.HTML(;
        canonical="https://cgt3.github.io/SearchableGeometries.jl",
        edit_link=get(ENV, "GITHUB_REF_TYPE", "") == "tag" ? :commit : "dev",
        assets=String[],
    ),
    pages=[
        "SearchableGeometries Documentation" => "index.md",
        "Geometry Types" => [
            "BoundingVolume" => "bounding_volume.md",
            "Ball" => "ball.md",
            "Hyperplane" => "hyperplane.md"
        ],
        "API Reference" => "api.md",
    ],
)

deploydocs(;
    repo="github.com/cgt3/SearchableGeometries.jl",
    devbranch="dev",
    devurl="dev"
)
