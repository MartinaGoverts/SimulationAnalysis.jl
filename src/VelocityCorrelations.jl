"""
    find_static_velocity_correlations(s::Simulation; kmin=7.0, kmax=7.4, kfactor=1)

Calculates the static velocity correlation function: ω(k) = 1/N < j(k)* j(k) >, 
where j(k) is a single-component current mode, the star represents 
complex conjugation and N is the total number of particles.

This function first constructs a `KSpace` object with the bounds specified by `kmin` and `kmax`,
then computes the current modes for this `KSpace`, and finally calculates the velocity correlation averaged over 
the k-vectors in the specified range.

# Arguments
- `s::Simulation`: A single- or multi-component simulation object.
- `kmin::Float64`: The lower bound over which to average ω(k).
- `kmax::Float64`: The upper bound over which to average ω(k).
- `kfactor = 1`: A multiplication factor for the construction of k-space.
"""
function find_static_velocity_correlations(s::Simulation; kmin=7.0, kmax=7.4, kfactor=1)
    kspace = construct_k_space(s, (kmin, kmax); kfactor=kfactor, negative=true, rectangular=false)
    jkt = find_current_modes(s, kspace; verbose=false)
    wk = find_static_velocity_correlations(s, kspace, jkt; kmin=kmin, kmax=kmax)
    return wk
end

"""
    find_static_velocity_correlations(s::Simulation, kspace::KSpace, jkt::AbstractDensityModes, k_sample_array::AbstractVector; k_binwidth=0.1)

Calculates the static velocity correlation function: ω(k) = 1/N < j(k)* j(k) >, with a pre-defined k-space
and current modes. ω(k) is calculated for each k-length in `k_sample_array`, and averaged over all k-vectors
within an interval `k_binwidth`.

# Arguments
- `s::Simulation`: A single- or multi-component simulation object.
- `kspace::KSpace`: A pre-computed `KSpace` object, which contains the available wavevectors.
- `jkt::AbstractDensityModes`: The real and imaginary parts of the pre-computed current modes.
- `k_sample_array::AbstractVector`: Contains the k-lengths at which to evaluate ω(k).
- `k_binwidth=0.1`: Specifies range of k-lengths over which to average (for every value in `k_sample_array`) 
"""
function find_static_velocity_correlations(s::Simulation, kspace::KSpace, jkt::AbstractDensityModes, k_sample_array::AbstractVector; k_binwidth=0.1)
    wk_array = []
    for (ik, k) in enumerate(k_sample_array)
        kmin = k - k_binwidth/2
        kmax = k + k_binwidth/2
        push!(wk_array, find_static_velocity_correlations(s, kspace, jkt; kmin=kmin, kmax=kmax))
    end
    return wk_array
end

"""
    find_static_velocity_correlations(s::Union{SingleComponentSimulation, SelfPropelledVoronoiSimulation}, kspace::KSpace, jkt::SingleComponentCurrentModes; kmin=0.0, kmax=10.0^10.0)

Calculates the static velocity correlation function for a single-component simulation. The velocity correlation function
is defined as: ω(k) = 1/N < j(k)* j(k) >.

ω(k) is calculated with pre-defined wavevectors and current modes, and averaged over all k-lengths between `kmin` and `kmax`.

# Arguments 
- `s::Union{SingleComponentSimulation, SelfPropelledVoronoiSimulation}`: A single-component simulation object.
- `kspace::KSpace`: A pre-computed `KSpace` object, which contains the available wavevectors.
- `jkt::SingleComponentCurrentModes`: The real and imaginary parts of the pre-computed single-component current modes.
- `kmin = 0.0`: The lower bound over which to average ω(k).
- `kmax = 10.0^10.0`: The upper bound over which to average ω(k).
"""
function find_static_velocity_correlations(s::Union{SingleComponentSimulation, SelfPropelledVoronoiSimulation}, kspace::KSpace, jkt::SingleComponentCurrentModes; kmin=0.0, kmax=10.0^10.0)
    wk = real_static_correlation_function(jkt.Re, jkt.Im, jkt.Re, jkt.Im, kspace, kmin, kmax)
    return wk / s.N
end

"""
    find_static_velocity_correlations(s::Union{MultiComponentSimulation,MCSPVSimulation}, kspace::KSpace, jkt::MultiComponentCurrentModes; kmin=0.0, kmax=10.0^10.0)

Calculates the partial static velocity correlation functions for a multicomponent simulation with pre-defined k-vectors and current modes, 
averaged over all k-vectors between `kmin` and `kmax`.
    
The partial velocity correlations are defined as: ω{αβ}(k) = 1/N < j{α}(k)* j{β}(k) >, where {.} denotes the particle species for which to calculate the correlation
function.

# Arguments
- `s::Union{MultiComponentSimulation,MCSPVSimulation}`: A multi-component simulation object.
- `kspace::KSpace`: A `KSpace` object, containing all available k-vectors.
- `jkt::MultiComponentCurrentModes`: The real and imaginary parts of the multi-component current modes.
- `kmin = 0.0`: The lower bound over which to average ω{αβ}(k).
- `kmax = 10.0^10.0`: The upper bound over which to average ω{αβ}(k).
"""
function find_static_velocity_correlations(s::Union{MultiComponentSimulation,MCSPVSimulation}, kspace::KSpace, jkt::MultiComponentCurrentModes; kmin=0.0, kmax=10.0^10.0)
    N_species = s.N_species
    wk = zeros(N_species, N_species)
    for α=1:N_species
        for β = 1:N_species
            wk[α,β] = real_static_correlation_function(jkt.Re[α], jkt.Im[α], jkt.Re[β], jkt.Im[β], kspace, kmin, kmax)
        end
    end
    return wk / s.N
end


"""
    find_dynamic_velocity_correlations(s::Simulation; kmin=7.0, kmax=7.4, kfactor=1)

Calculates the dynamic velocity correlation function  ω(k,t) = 1/N < j(k)* j(k,t) >,
where j(k,t) is a single-component current mode, the star represents complex conjugation and N is the total number of particles.

This function first constructs a `KSpace` object with the bounds specified by `kmin` and `kmax`,
then computes the current modes for this `KSpace`, and finally calculates the velocity correlation averaged over 
the k-vectors in the specified range. ω(k,t) is calculated for all time intervals in the simulation `s`. Note that ω(k,t) typically decays quickly, so the data should contain small time intervals to obtain satisfactory results.

# Arguments
- `s::Simulation`: A simulation object.
- `kmin=7.0`: The lower bound from which to calculate the `KSpace` object.
- `kmax=7.4`: The upper bound from which to calculate the `KSpace` object.
- `kfactor=1`: A multiplication factor for the construction of k-space.
"""
function find_dynamic_velocity_correlations(s::Simulation; kmin=7.0, kmax=7.4, kfactor=1)
    kspace = construct_k_space(s, (kmin, kmax); kfactor=kfactor, negative=true, rectangular=false)
    jkt = find_current_modes(s, kspace; verbose=false)
    wkt = find_dynamic_velocity_correlations(s, kspace, jkt; kmin=kmin, kmax=kmax)
    return wkt
end

"""
    find_dynamic_velocity_correlations(s::Simulation, kspace::KSpace, jkt::AbstractDensityModes, k_sample_array::AbstractVector; k_binwidth=0.1)

Calculates the dynamic velocity correlation function  ω(k,t) = 1/N < j(k)* j(k,t) >, for a pre-defined space of k-vectors and current modes.
ω(k,t) is calculated for every k-length in `k_sample_array`, and averaged over k-vectors in a range `k_binwidth` of each value in `k_sample_array`.
The dynamic velocity correlations are calculated for every time interval in `s`.

# Arguments
- `s::Simulation`: A simulation object.
- `kspace::KSpace`: Contains the available wavevectors.
- `jkt::AbstractDensityModes`: The single- or multi-component current modes.
- `k_sample_array`: Contains the wavevectors at which to evaluate ω(k,t).
- `kmin=7.0`: The lower bound over which to average ω(k,t).
- `kmax=7.4`: The upper bound over which to average ω(k,t).
"""
function find_dynamic_velocity_correlations(s::Simulation, kspace::KSpace, jkt::AbstractDensityModes, k_sample_array::AbstractVector; k_binwidth=0.1)
    wkt_array = []
    for (ik, k) in enumerate(k_sample_array)
        kmin = k - k_binwidth/2
        kmax = k + k_binwidth/2
        push!(wkt_array, find_dynamic_velocity_correlations(s, kspace, jkt; kmin=kmin, kmax=kmax))
    end
    return wkt_array
end

"""
    find_dynamic_velocity_correlations(s::Union{SingleComponentSimulation, SelfPropelledVoronoiSimulation}, kspace::KSpace, jkt::SingleComponentCurrentModes; kmin=0.0, kmax=10.0^10.0)

Calculates the dynamic velocity correlations for a single-component simulation, for all time intervals and with a predefined kspace and current modes.
The dynamic velocity correlation function is defined as ω(k,t) = 1/N < j(k)* j(k,t) >. The result is averaged over all wavevectors between `kmin` and `kmax`.

# Arguments
- `s::Union{SingleComponentSimulation, SelfPropelledVoronoiSimulation}`: A single-component simulation object.
- `kspace::KSpace`: Contains the available wavevectors.
- `jkt::SingleComponentCurrentModes`: The single-component current modes.
- `kmin=7.0`: The lower bound over which to average ω(k,t).
- `kmax=7.4`: The upper bound over which to average ω(k,t).
"""
function find_dynamic_velocity_correlations(s::Union{SingleComponentSimulation, SelfPropelledVoronoiSimulation}, kspace::KSpace, jkt::SingleComponentCurrentModes; kmin=0.0, kmax=10.0^10.0)
    Ndt = length(s.dt_array)
    wkt = zeros(Ndt)
    real_correlation_function!(wkt, jkt.Re, jkt.Im, jkt.Re, jkt.Im, kspace, s.dt_array, s.t1_t2_pair_array, kmin, kmax)
    return wkt ./ s.N
end

"""
    find_dynamic_velocity_correlations(s::Union{MultiComponentSimulation,MCSPVSimulation}, kspace::KSpace, jkt::MultiComponentCurrentModes; kmin=0.0, kmax=10.0^10.0)

Calculates the partial dynamic velocity correlations for all species and for all time intervals in a multi-component simulation.
The partial dynamic velocity correlation function is defined as ω{αβ}(k,t) = 1/N < j{α}(k)* j{β}(k,t) , where {.} denotes the species for which to evaluate.
The result is averaged over all k-vectors with magnitude between `kmin` and `kmax`.

# Arguments
- `s::Union{MultiComponentSimulation, MCSPVSimulation}`: A multi-component simulation object.
- `kspace::KSpace`: Contains the available wavevectors.
- `jkt::MultiComponentCurrentModes`: The multi-component current modes.
- `kmin=0.0`: The lower bound over which to average ω{αβ}(k,t).
- `kmax=10.0^10.0`: The upper bound over which to average ω{αβ}(k,t).
"""
function find_dynamic_velocity_correlations(s::Union{MultiComponentSimulation,MCSPVSimulation}, kspace::KSpace, jkt::MultiComponentCurrentModes; kmin=0.0, kmax=10.0^10.0)
    N_species = s.N_species
    Ndt = length(s.dt_array)
    wkt = [zeros(Ndt) for α=1:N_species, β=1:N_species]
    for α=1:N_species
        for β = 1:N_species
            real_correlation_function!(wkt[α,β], jkt.Re[α], jkt.Im[α], jkt.Re[β], jkt.Im[β], kspace, s.dt_array, s.t1_t2_pair_array, kmin, kmax)
        end
    end
    return wkt ./ s.N
end
