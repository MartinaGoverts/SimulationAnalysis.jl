# SimulationAnalysis.jl

[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://IlianPihlajamaa.github.io/SimulationAnalysis.jl/dev)

`SimulationAnalysis.jl` is a Julia package for analyzing data from particle-based simulations. It provides a collection of tools for computing various physical quantities and performing common analysis tasks.

## Features

*   Load simulation data from some limited formats.
*   Calculate mean squared displacement.
*   Compute radial distribution functions.
*   Analyze structure factors and intermediate scattering functions.
*   Perform clustering analysis.

## Installation

The package can be installed from the Julia package manager. From the Julia REPL, type `]` to enter Pkg mode, and then run:

```julia
pkg> add "https://github.com/IlianPihlajamaa/SimulationAnalysis.jl"
```

## Usage

Here is a basic example of how to use the package to load a simulation and calculate the mean squared displacement:

```julia
using SimulationAnalysis

# Load simulation data
file = joinpath("test", "data", "test_trajectory.h5")
simulation = SimulationAnalysis.read_continuously_hard_sphere_simulation(file; time_origins=10)

# Calculate the mean squared displacement
msd = find_mean_squared_displacement(simulation)

println("Mean Squared Displacement:")
println(msd)
```

For more detailed information and a complete list of features, please see the [documentation](https://IlianPihlajamaa.github.io/SimulationAnalysis.jl/dev).
