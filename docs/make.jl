using Documenter
using SimulationAnalysis

makedocs(
    sitename = "SimulationAnalysis.jl",
    format = Documenter.HTML(),
    modules = [SimulationAnalysis],
    pages = [
        "Home" => "index.md",
        "Tutorial" => "tutorial.md",
        "Examples" => "examples.md",
        "API" => "api.md",
    ],
    doctest = false
)

deploydocs(
    repo = "github.com/IlianPihlajamaa/SimulationAnalysis.jl.git",
)
