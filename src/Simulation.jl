
"""
    SingleComponentSimulation

A struct to hold the data of a single-component simulation.

# Fields
- `N::Int64`: Number of particles.
- `Ndims::Int64`: Number of dimensions.
- `r_array::Array{Float64, 3}`: Particle positions.
- `v_array::Array{Float64, 3}`: Particle velocities.
- `F_array::Array{Float64, 3}`: Particle forces.
- `D_array::Array{Float64, 1}`: Diffusion coefficient.
- `t_array::Array{Float64, 1}`: Time steps.
- `box_sizes::Vector{Float64}`: Box sizes.
- `dt_array::Array{Int64, 1}`: Time step differences.
- `t1_t2_pair_array::Vector{Array{Int64, 2}}`: Time pairs.
- `filepath::String`: Path to the simulation file.
"""
mutable struct SingleComponentSimulation <: Simulation
    N::Int64
    Ndims::Int64
    r_array::Array{Float64, 3}
    v_array::Array{Float64, 3}
    F_array::Array{Float64, 3}
    D_array::Array{Float64, 1}
    t_array::Array{Float64, 1}
    box_sizes::Vector{Float64}
    dt_array::Array{Int64, 1}
    t1_t2_pair_array::Vector{Array{Int64, 2}}
    filepath::String
end

"""
    MultiComponentSimulation

A struct to hold the data of a multi-component simulation.

# Fields
- `N::Int64`: Total number of particles.
- `Ndims::Int64`: Number of dimensions.
- `N_species::Int64`: Number of species.
- `N_particles_per_species::Vector{Int}`: Number of particles per species.
- `r_array::Vector{Array{Float64, 3}}`: Particle positions for each species.
- `v_array::Vector{Array{Float64, 3}}`: Particle velocities for each species.
- `F_array::Vector{Array{Float64, 3}}`: Particle forces for each species.
- `t_array::Vector{Float64}`: Time steps.
- `box_sizes::Vector{Float64}`: Box sizes.
- `dt_array::Vector{Int}`: Time step differences.
- `t1_t2_pair_array::Vector{Array{Int64, 2}}`: Time pairs.
- `filepath::String`: Path to the simulation file.
"""
mutable struct MultiComponentSimulation <: Simulation
    N::Int64
    Ndims::Int64
    N_species::Int64
    N_particles_per_species::Vector{Int}
    r_array::Vector{Array{Float64, 3}}
    v_array::Vector{Array{Float64, 3}}
    F_array::Vector{Array{Float64, 3}}
    t_array::Vector{Float64}
    box_sizes::Vector{Float64}
    dt_array::Vector{Int}
    t1_t2_pair_array::Vector{Array{Int64, 2}}
    filepath::String
end