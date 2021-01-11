#julia -e 'using LiveServer; serve(dir="./build")'
#julia make.jl

using Documenter, ProbabilityBoundsAnalysis

makedocs(
    modules = [ProbabilityBoundsAnalysis],
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    sitename="ProbabilityBoundsAnalysis.jl",
    authors = "Ander Gray and Scott Ferson",
    doctest=false,
    pages = [
        "Home" => "index.md",
        "Theory" => "theory.md",
        "Constructing p-boxes" => "pboxes.md",
        "Arithmetic" => "arithmetic.md",
        "Dependence" => "copulas.md",
        "User guide" => "userguide.md",
        "Examples" => "examples.md",
        "API" => "api.md"
        ]
    )



deploydocs(
    repo   = "github.com/AnderGray/ProbabilityBoundsAnalysis.jl",
)