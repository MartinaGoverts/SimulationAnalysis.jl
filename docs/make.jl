using Documenter
using SimulationAnalysis

makedocs(
    sitename = "SimulationAnalysis.jl",
    format = Documenter.HTML(),
    modules = [SimulationAnalysis],
    pages = [
        "Home" => "index.md",
        "API" => "api.md",
        "Tutorial" => "tutorial.md",
        "Examples" => "examples.md"
    ],
    doctest = false
)

deploydocs(
    repo = "github.com/IlianPihlajamaa/SimulationAnalysis.jl.git",
)
