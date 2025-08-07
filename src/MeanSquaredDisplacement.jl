"""
    find_mean_squared_displacement(simulation::Simulation; per_particle=false)

Calculates the mean squared displacement for a simulation.

# Arguments
- `simulation::Simulation`: The simulation.
- `per_particle=false`: Whether to return the mean squared displacement per particle.

# Returns
- `msd`: The mean squared displacement.
- `msd_per_particle`: The mean squared displacement per particle (if `per_particle` is true).
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
    find_mean_squared_displacement(r, dt_array, t1_t2_pair_array; verbose=true)

Calculates the mean squared displacement from the t1, t2 pairs.

# Arguments
- `r`: The positions of the particles.
- `dt_array`: The array of time differences.
- `t1_t2_pair_array`: The array of pairs of times.
- `verbose=true`: Whether to print verbose output.

# Returns
- `msd`: The mean squared displacement.
- `msd_per_particle`: The mean squared displacement per particle.
"""
function find_mean_squared_displacement(r, dt_array, t1_t2_pair_array; verbose=true)
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
