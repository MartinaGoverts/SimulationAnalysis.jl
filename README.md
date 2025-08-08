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

Here is a basic example of how to use the package to load a simulation and calculate radial distribution function:

```julia
using  SimulationAnalysis
# Set parameters for g(r) calculation
Nbins = 100
rmax = 5.0

# Calculate the radial distribution function
r, g_r = SimulationAnalysis.find_radial_distribution_function(sim, Nbins, rmax)

using Plots
# Plot the g(r)
plot(r, g_r,
    xlabel="r",
    ylabel="g(r)",
    title="Radial Distribution Function",
    legend=false,
    lw=2
)
```

![gr](docs/src/plots/gr.png)

For more detailed information and a complete list of features, please see the [documentation](https://IlianPihlajamaa.github.io/SimulationAnalysis.jl/dev).
