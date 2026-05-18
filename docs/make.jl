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
    pagesonly=true,
    checkdocs=:none,
    format=Documenter.HTML(;
        canonical="https://cgt3.github.io/SearchableGeometries.jl",
        edit_link=get(ENV, "GITHUB_REF_TYPE", "") == "tag" ? :commit : "dev",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Manual" => [
            "Getting Started" => "man/getting_started.md",
            "Geometric Search" => "man/geometric_search.md",
            "Bounding Volumes" => "man/bounding_volumes.md",
            "Balls" => "man/balls.md",
            "Hyperplanes" => "man/hyperplanes.md",
            "Cones" => "man/cones.md",
        ], "Library" => [
            "API" => "lib/api.md",
            "Geometry Types" => [
                "Bounding Volumes" => "lib/types/bounding_volumes.md",
                "Balls" => "lib/types/balls.md",
                "Hyperplanes" => "lib/types/hyperplanes.md",
            ],
            # "Operations" => [
            #     "Containment" => "lib/operations/containment.md",
            #     "Intersection" => "lib/operations/intersection.md",
            #     "Closest and Furthest Points" => "lib/operations/closest_and_furthest_points.md"
            # ],
        ],
    ],
)

deploydocs(;
    repo="github.com/cgt3/SearchableGeometries.jl",
    devbranch="dev",
    devurl="dev"
)
