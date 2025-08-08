
"""
    SingleComponentSimulation

A mutable struct that holds all data for a single-component particle simulation.

# Fields
- `N::Int64`: Number of particles.
- `Ndims::Int64`: Number of spatial dimensions.
- `Nt::Int64`: Number of time steps.
- `dt::Float64`: Time step size.
- `r_array::Array{Float64, 3}`: Particle positions. The dimensions of the array are `(N, Ndims, N_time_steps)`.
- `v_array::Array{Float64, 3}`: Particle velocities. The dimensions of the array are `(N, Ndims, N_time_steps)`.
- `F_array::Array{Float64, 3}`: Forces on particles. The dimensions of the array are `(N, Ndims, N_time_steps)`.
- `D_array::Array{Float64, 1}`: Diameter for each particle.
- `t_array::Array{Float64, 1}`: Vector of time points corresponding to each time step.
- `box_sizes::Vector{Float64}`: A vector holding the simulation box dimensions, e.g., `[Lx, Ly, Lz]`.
- `dt_array::Array{Float64, 1}`: An array of time step differences (`Δt`) used for calculating time correlation functions.
- `t1_t2_pair_array::Vector{Array{Int64, 2}}`: An array of `(t1, t2)` index pairs, used for efficient calculation of time correlations.
- `filepath::String`: Path to the file from which the simulation data was loaded.
"""
mutable struct SingleComponentSimulation <: Simulation
    N::Int64
    Ndims::Int64
    Nt::Int64
    dt::Float64
    r_array::Array{Float64, 3}
    v_array::Array{Float64, 3}
    F_array::Array{Float64, 3}
    D_array::Vector{Float64}
    t_array::Vector{Float64}
    box_sizes::Vector{Float64}
    dt_array::Vector{Float64}
    t1_t2_pair_array::Vector{Array{Int64, 2}}
    filepath::String
end

"""
    MultiComponentSimulation

A mutable struct that holds all data for a multi-component particle simulation.

# Fields
- `N::Int64`: Total number of particles across all species.
- `Ndims::Int64`: Number of spatial dimensions.
- `N_species::Int64`: Number of different particle species.
- `Nt::Int64`: Number of time steps.
- `dt::Float64`: Time step size.
- `N_particles_per_species::Vector{Int}`: A vector containing the number of particles for each species.
- `r_array::Vector{Array{Float64, 3}}`: A vector of particle position arrays. Each element `r_array[i]` is a `(N_particles_per_species[i], Ndims, N_time_steps)` array for species `i`.
- `v_array::Vector{Array{Float64, 3}}`: A vector of particle velocity arrays for each species.
- `F_array::Vector{Array{Float64, 3}}`: A vector of particle force arrays for each species.
- `t_array::Vector{Float64}`: Vector of time points corresponding to each time step.
- `box_sizes::Vector{Float64}`: A vector holding the simulation box dimensions, e.g., `[Lx, Ly, Lz]`.
- `dt_array::Vector{Float64}`: An array of time step differences (`Δt`) used for calculating time correlation functions.
- `t1_t2_pair_array::Vector{Array{Int64, 2}}`: An array of `(t1, t2)` index pairs, used for efficient calculation of time correlations.
- `filepath::String`: Path to the file from which the simulation data was loaded.
"""
mutable struct MultiComponentSimulation <: Simulation
    N::Int64
    Ndims::Int64
    N_species::Int64
    Nt::Int64
    dt::Float64
    N_particles_per_species::Vector{Int}
    r_array::Vector{Array{Float64, 3}}
    v_array::Vector{Array{Float64, 3}}
    F_array::Vector{Array{Float64, 3}}
    t_array::Vector{Float64}
    box_sizes::Vector{Float64}
    dt_array::Vector{Float64}
    t1_t2_pair_array::Vector{Array{Int64, 2}}
    filepath::String
end

"""
    SelfPropelledVoronoiSimulation

A mutable struct that holds all data for a self-propelled Voronoi particle simulation.

# Fields
- `N::Int64`: Number of particles.
- `Ndims::Int64`: Number of spatial dimensions.
- `Nt::Int64`: Number of time steps.
- `dt::Float64`: Time step size.
- `r_array::Array{Float64, 3}`: Particle positions.
- `u_array::Array{Float64, 2}`: Particle orientations (radians).
- `F_array::Array{Float64, 3}`: Forces on particles.
- `perimeter_array::Array{Float64, 2}`: Perimeters of Voronoi cells.
- `area_array::Array{Float64, 2}`: Areas of Voronoi cells.
- `Epot_array::Array{Float64, 1}`: Potential energies.
- `t_array::Vector{Float64}`: Vector of time points.
- `box_sizes::Vector{Float64}`: Dimensions of the simulation box.
- `dt_array::Vector{Float64}`: Time step sizes.
- `t1_t2_pair_array::Vector{Array{Int64, 2}}`: Array of `(t1, t2)` index pairs.
- `filepath::String`: Path to the file from which the simulation data was loaded.
"""
struct SelfPropelledVoronoiSimulation <: Simulation
    N::Int64
    Ndims::Int64
    Nt::Int64
    dt::Float64
    r_array::Array{Float64, 3}
    u_array::Array{Float64, 2}
    F_array::Array{Float64, 3}
    perimeter_array::Array{Float64, 2}
    area_array::Array{Float64, 2}
    Epot_array::Vector{Float64} 
    t_array::Vector{Float64}
    box_sizes::Vector{Float64}
    dt_array::Vector{Float64}
    t1_t2_pair_array::Vector{Array{Int64, 2}}
    filepath::String
end