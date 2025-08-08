"""
    find_mean_squared_displacement(simulation::Simulation; per_particle=false)

Calculates the mean squared displacement (MSD) for a given simulation.

The MSD is defined as `<|r(t) - r(0)|^2>`, where the average is taken over all particles and time origins.

# Arguments
- `simulation::Simulation`: The simulation data.
- `per_particle::Bool=false`: If `true`, the function returns the MSD calculated for each particle individually, in addition to the total MSD.
- `power::Int=2`: The power to which the displacement is raised. Default is 2 for MSD. Set to 4 for the mean quartic displacement of the displacement for example.

# Returns
If `per_particle` is `false` (default):
- `msd::Vector{Float64}`: A vector containing the MSD for each time delay `Δt` in `simulation.dt_array`.

If `per_particle` is `true`:
- `msd::Vector{Float64}`: The total MSD vector.
- `msd_per_particle::Matrix{Float64}`: A `(N, N_dt)` matrix of MSDs for each particle.
"""
function find_mean_squared_displacement(simulation::Simulation; per_particle=false, power=2)
    msd, msd_per_particle = find_mean_squared_displacement(simulation.r_array, simulation.dt_array, simulation.t1_t2_pair_array, power)
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
- `power::Int=2`: The power to which the displacement is raised. 

# Returns
- `msd::Vector{Float64}`: A vector containing the MSD for each time delay `Δt` in `dt_array`.
- `msd_per_particle::Matrix{Float64}`: A `(N, N_dt)` matrix of MSDs for each particle.
"""
function find_mean_squared_displacement(r, dt_array, t1_t2_pair_array, power)
    dims , N, _ = size(r)  
    N_dt = length(dt_array)
    msd_per_particle = zeros(N, N_dt)
    pd2 = power / 2
    for particle = 1:N
        for idt = eachindex(dt_array)
            msddtp1 = 0.0
            pairs_idt = t1_t2_pair_array[idt]
            Npairs = size(pairs_idt,1)
            for ipair = 1:Npairs
                t1 = pairs_idt[ipair,1]
                t2 = pairs_idt[ipair,2]
                r2abs = 0.0
                for dim = 1:dims
                    r2abs += (r[dim, particle, t2] - r[dim, particle, t1])^2
                end
                if power == 2
                    msddtp1 += r2abs
                else
                    msddtp1 += r2abs ^ pd2
                end
            end
            msd_per_particle[particle, idt] = msddtp1 / Npairs
        end
    end
    msd = reshape(sum(msd_per_particle, dims=1), N_dt) 
    msd .= msd ./ N

    return msd, msd_per_particle
end 

"""
    α2(msd, mqd, dims)
Calculates the non-Gaussian parameter α2 from the mean squared displacement (MSD) and mean quartic displacement (MQD).
This parameter is used to quantify the deviation from Gaussian behavior in particle displacements.
# Arguments
- `msd::Float64`: The mean squared displacement.
- `mqd::Float64`: The mean quartic displacement.
- `dims::Int`: The number of spatial dimensions in the simulation.
# Returns
- `α2::Float64`: The non-Gaussian parameter.
"""
α2(msd, mqd, dims) = dims/(dims+2) * mqd / msd^2 - 1

"""
    find_non_gaussian_parameter(simulation::Simulation; per_particle=false)
"""
function find_non_gaussian_parameter(simulation::Simulation; per_particle=false)
    msd, msd_per_particle = find_mean_squared_displacement(simulation.r_array, simulation.dt_array, simulation.t1_t2_pair_array, 2)
    mqd, mqd_per_particle = find_mean_squared_displacement(simulation.r_array, simulation.dt_array, simulation.t1_t2_pair_array, 4)
    dims = simulation.Ndims
    alpha2 = α2.(msd, mqd, dims)
    if per_particle
        alpha2_per_particle = α2.(msd_per_particle, mqd_per_particle, dims)
        return alpha2, alpha2_per_particle
    else
        return alpha2
    end
end

