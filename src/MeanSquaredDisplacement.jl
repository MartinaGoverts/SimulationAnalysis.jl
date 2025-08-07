"""
    find_mean_squared_displacement(simulation::Simulation; per_particle=false)

Calculates the mean squared displacement (MSD) for a given simulation.

The MSD is defined as `<|r(t) - r(0)|^2>`, where the average is taken over all particles and time origins.

# Arguments
- `simulation::Simulation`: The simulation data.
- `per_particle::Bool=false`: If `true`, the function returns the MSD calculated for each particle individually, in addition to the total MSD.

# Returns
If `per_particle` is `false` (default):
- `msd::Vector{Float64}`: A vector containing the MSD for each time delay `Δt` in `simulation.dt_array`.

If `per_particle` is `true`:
- `msd::Vector{Float64}`: The total MSD vector.
- `msd_per_particle::Matrix{Float64}`: A `(N, N_dt)` matrix of MSDs for each particle.
"""
function find_mean_squared_displacement(simulation::Simulation; per_particle=false)
    msd, msd_per_particle = find_mean_squared_displacement(simulation.r_array, simulation.dt_array, simulation.t1_t2_pair_array)
    if per_particle
        return msd, msd_per_particle
    else
        return msd
    end
end


"""
    find_mean_squared_displacement(r, dt_array, t1_t2_pair_array)

Calculates the mean squared displacement (MSD) using pre-calculated time pairs.

This is the core implementation for calculating the MSD. It iterates through particles and time pairs to compute the displacement.

# Arguments
- `r::Array{Float64, 3}`: A `(Ndims, N, N_timesteps)` array of particle positions.
- `dt_array::Vector{Int}`: An array of time step differences (`Δt`).
- `t1_t2_pair_array::Vector{Array{Int64, 2}}`: An array of `(t1, t2)` index pairs for each `Δt`.

# Returns
- `msd::Vector{Float64}`: A vector containing the MSD for each time delay `Δt` in `dt_array`.
- `msd_per_particle::Matrix{Float64}`: A `(N, N_dt)` matrix of MSDs for each particle.
"""
function find_mean_squared_displacement(r, dt_array, t1_t2_pair_array)
    dims , N, _ = size(r)  
    N_dt = length(dt_array)
    msd_per_particle = zeros(N, N_dt)
    for particle = 1:N
        for idt = eachindex(dt_array)
            msddtp1 = 0.0
            pairs_idt = t1_t2_pair_array[idt]
            Npairs = size(pairs_idt,1)
            for ipair = 1:Npairs
                t1 = pairs_idt[ipair,1]
                t2 = pairs_idt[ipair,2]
                for dim = 1:dims
                    msddtp1 += (r[dim, particle, t2] - r[dim, particle, t1])^2
                end
            end
            msd_per_particle[particle, idt] = msddtp1
        end
    end
    msd = reshape(sum(msd_per_particle, dims=1), N_dt)
    for idt in eachindex(dt_array)
        pairs_idt = t1_t2_pair_array[idt]
        Npairs = size(pairs_idt,1)
        msd[idt] /= N*Npairs
        msd_per_particle[:, idt] ./= Npairs
    end
    return msd, msd_per_particle
end 
