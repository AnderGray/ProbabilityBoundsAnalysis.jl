#julia -e 'using LiveServer; serve(dir="./build")'
#julia make.jl

using ProbabilityBoundsAnalysis, IntervalArithmetic, Documenter

DocMeta.setdocmeta!(ProbabilityBoundsAnalysis, :DocTestSetup, :(using ProbabilityBoundsAnalysis, IntervalArithmetic, Random; Random.seed!(MersenneTwister(), 31415)); recursive=true)

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
        "Set operations and comparisons" => "sets.md",
        "Examples" => "examples.md",
        "API" => "api.md"
        ]
    )



deploydocs(
    repo   = "github.com/AnderGray/ProbabilityBoundsAnalysis.jl",
)
