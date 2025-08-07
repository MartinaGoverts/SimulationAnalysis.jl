# Examples

Here are some examples of how to use `SimulationAnalysis.jl`.

## Mean Squared Displacement

This example shows how to calculate the mean squared displacement (MSD) for a simulation.

```julia
using SimulationAnalysis

# Load a simulation
filename = "path/to/your/simulation.dat"
dt = 0.005
sim = read_WCA_simulation(filename, dt)

# Calculate the MSD
msd = find_mean_squared_displacement(sim)

# The `msd` variable now contains the mean squared displacement as a function of time.
```

## Intermediate Scattering Function

This example shows how to calculate the intermediate scattering function (ISF) for a simulation.

```julia
using SimulationAnalysis

# Load a simulation
filename = "path/to/your/simulation.dat"
dt = 0.005
sim = read_WCA_simulation(filename, dt)

# Define the k-range for the ISF
kmin = 7.0
kmax = 7.4

# Calculate the ISF
isf = find_intermediate_scattering_function(sim; kmin=kmin, kmax=kmax)

# The `isf` variable now contains the intermediate scattering function as a function of time.
```
