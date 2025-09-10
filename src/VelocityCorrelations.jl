# TO DO: dynamic velocity correlations (like the intermediate scattering function calculation)

"""
    find_static_velocity_correlations(s::Simulation; kmin=7.0, kmax=7.4, kfactor=1)

[...]
"""
function find_static_velocity_correlations(s::Simulation; kmin=7.0, kmax=7.4, kfactor=1)
    kspace = construct_k_space(s, (kmin, kmax); kfactor=kfactor, negative=true, rectangular=false)
    jkt = find_current_modes(s, kspace; verbose=false)
    wk = find_static_velocity_correlations(s, kspace, jkt; kmin=kmin, kmax=kmax)
    return wk
end

"""
    find_structure_factor(s::Simulation, kspace::KSpace, ρkt::AbstractDensityModes, k_sample_array::AbstractVector; k_binwidth=0.1)

[...]
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
    find_structure_factor(s::Union{SingleComponentSimulation, SelfPropelledVoronoiSimulation}, kspace::KSpace, jkt::SingleComponentCurrentModes; kmin=0.0, kmax=10.0^10.0)

[...]
"""
function find_static_velocity_correlations(s::Union{SingleComponentSimulation, SelfPropelledVoronoiSimulation}, kspace::KSpace, jkt::SingleComponentCurrentModes; kmin=0.0, kmax=10.0^10.0)
    wk = real_static_correlation_function(jkt.Re, jkt.Im, jkt.Re, jkt.Im, kspace, kmin, kmax)
    return wk / s.N
end

"""
    find_structure_factor(s::Union{MultiComponentSimulation,MCSPVSimulation}, kspace::KSpace, jkt::MultiComponentCurrentModes; kmin=0.0, kmax=10.0^10.0)

[...]
"""
function find_static_velocity_correlations(s::Union{MultiComponentSimulation,MCSPVSimulation}, kspace::KSpace, jkt::MultiComponentCurrentModes; kmin=0.0, kmax=10.0^10.0)
    N_species = s.N_species
    wk = zeros(N_species, N_species)
    for α=1:N_species
        for β = α:N_species  # order of alpha, beta?
            wk[β,α] = real_static_correlation_function(jkt.Re[α], jkt.Im[α], jkt.Re[β], jkt.Im[β], kspace, kmin, kmax)
            if α != β
                wk[α, β] = wk[β, α]
            end
        end
    end
    return wk / s.N
end

