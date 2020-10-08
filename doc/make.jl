#julia -e 'using LiveServer; serve(dir="./build")'
#julia make.jl

using Documenter, ProbabilityBoundsAnalysis

makedocs(
    modules = [ProbabilityBoundsAnalysis],
    sitename="ProbabilityBoundsAnalysis.jl",
    authors = "Ander Gray",
    pages = ["Home" => "index.md"]
    )