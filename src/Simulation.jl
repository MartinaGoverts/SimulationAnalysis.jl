
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