# API Reference

## Simulation

```@docs
SimulationAnalysis.SingleComponentSimulation
SimulationAnalysis.MultiComponentSimulation
```

## Load Data

```@docs
SimulationAnalysis.read_WCA_simulation
SimulationAnalysis.read_Newtonian_KAWCA_simulation
SimulationAnalysis.read_Brownian_KALJ_simulation
SimulationAnalysis.read_Brownian_KAWCA_simulation
SimulationAnalysis.read_monodisperse_hard_sphere_simulation
SimulationAnalysis.read_continuously_hard_sphere_simulation
```

## K-Space

```@docs
SimulationAnalysis.KSpace
SimulationAnalysis.construct_k_space
```

## Correlation Function

```@docs
SimulationAnalysis.find_correlation_function
SimulationAnalysis.find_correlation_matrix
```

## Density Modes

```@docs
SimulationAnalysis.find_density_modes
SimulationAnalysis.SingleComponentDensityModes
SimulationAnalysis.MultiComponentDensityModes
```

## F4 Diagonal

```@docs
SimulationAnalysis.find_F4_diagonal
SimulationAnalysis.find_F4_diagonal_all_k
SimulationAnalysis.find_F4_super_diagonal
```

## Forces

```@docs
SimulationAnalysis.Weysser
SimulationAnalysis.Berthier
SimulationAnalysis.WCA
SimulationAnalysis.KAWCA
SimulationAnalysis.calculate_forces!
```

## Intermediate Scattering Function

```@docs
SimulationAnalysis.find_intermediate_scattering_function
SimulationAnalysis.find_self_intermediate_scattering_function
```

## Structure Factors

```@docs
SimulationAnalysis.find_structure_factor
SimulationAnalysis.find_S4_offiagonal
```

## Mean Squared Displacement

```@docs
SimulationAnalysis.find_mean_squared_displacement
```

## Overlap Function

```@docs
SimulationAnalysis.find_overlap_function
```

## Neighborlists

```@docs
SimulationAnalysis.find_relative_distance_neighborlists
SimulationAnalysis.find_absolute_distance_neighborlists
SimulationAnalysis.find_voronoi_neighborlists
```

## Bond Breaking Parameter

```@docs
SimulationAnalysis.find_CB
SimulationAnalysis.find_χBB_smoothed
SimulationAnalysis.find_chi_BB
```

## Utils

```@docs
SimulationAnalysis.find_relaxation_time
```

## Radial Distribution Function

```@docs
SimulationAnalysis.find_radial_distribution_function
```
