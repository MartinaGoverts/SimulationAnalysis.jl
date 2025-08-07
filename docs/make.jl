using Documenter
using SimulationAnalysis

makedocs(
    sitename = "SimulationAnalysis.jl",
    format = Documenter.HTML(),
    modules = [SimulationAnalysis]
)

deploydocs(
    repo = "github.com/ipihlama/SimulationAnalysis.jl.git",
)
