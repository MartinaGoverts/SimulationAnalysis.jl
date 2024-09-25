function find_intermediate_scattering_function(s::Simulation; kmin=7.0, kmax=7.4, kfactor=1)
    kspace = construct_k_space(s, (kmin, kmax); kfactor=kfactor, negative=true, rectangular=false)
    ρkt = find_density_modes(s, kspace; verbose=false)
    F = find_intermediate_scattering_function(s, kspace, ρkt; kmin=kmin, kmax=kmax)
    return F
end

function find_intermediate_scattering_function(s::Simulation, kspace::KSpace, ρkt, k_sample_array::AbstractVector; k_binwidth=0.1)
    F_array = []
    for (ik, k) in enumerate(k_sample_array)
        kmin = k - k_binwidth/2
        kmax = k + k_binwidth/2
        push!(F_array, find_intermediate_scattering_function(s, kspace, ρkt; kmin=kmin, kmax=kmax))
    end
    return F_array
end

function find_intermediate_scattering_function(s::SingleComponentSimulation, kspace::KSpace, ρkt::SingleComponentDensityModes; kmin=0.0, kmax=10.0^10.0)
    Ndt = length(s.dt_array)
    Fk = zeros(Ndt)
    real_correlation_function!(Fk, ρkt.Re, ρkt.Im, ρkt.Re, ρkt.Im, kspace, s.dt_array, s.t1_t2_pair_array, kmin, kmax)
    return Fk ./ s.N
end

function find_intermediate_scattering_function(s::MultiComponentSimulation, kspace::KSpace, ρkt::MultiComponentDensityModes; kmin=0.0, kmax=10.0^10.0)
    N_species = s.N_species
    Ndt = length(s.dt_array)
    Fk = [zeros(Ndt) for α=1:N_species, β=1:N_species]
    for α=1:N_species
        for β = 1:N_species
            real_correlation_function!(Fk[α,β], ρkt.Re[α], ρkt.Im[α], ρkt.Re[β], ρkt.Im[β], kspace, s.dt_array, s.t1_t2_pair_array, kmin, kmax)
            # Fk[β, α] .= Fk[α,β]
        end
    end
    return Fk ./ s.N
end

function find_self_intermediate_scattering_function(s::Simulation; kmin=7.0, kmax=7.4, kfactor=1)
    kspace = construct_k_space(s, (kmin, kmax); kfactor=kfactor, negative=true, rectangular=false)
    F = find_self_intermediate_scattering_function(s, kspace; kmin=kmin, kmax=kmax)
    return F
end

function find_self_intermediate_scattering_function(s::Simulation, kspace::KSpace, k_sample_array::AbstractVector; k_binwidth=0.1)
    F_array = []
    for (ik, k) in enumerate(k_sample_array)
        kmin = k - k_binwidth/2
        kmax = k + k_binwidth/2
        push!(F_array, find_self_intermediate_scattering_function(s, kspace; kmin=kmin, kmax=kmax))
    end
    return F_array
end

function find_self_intermediate_scattering_function(s::SingleComponentSimulation, kspace::KSpace; kmin=0.0, kmax=10.0^10.0)
    kmask = (kmin .< kspace.k_lengths .< kmax)
    k_lengths = kspace.k_lengths[kmask]
    k_array = kspace.k_array[:, kmask]
    N = s.N
    Nk = length(k_lengths)
    Ndt = length(s.dt_array)
    Nt = length(s.t_array)
    coskr = zeros(Nt, Nk)
    sinkr = zeros(Nt, Nk)
    Fks_per_particle = zeros(Ndt, N)
    for particle = 1:N
        if s.Ndims == 3
            @turbo for ik = 1:Nk
                kx = k_array[1, ik]
                ky = k_array[2, ik]
                kz = k_array[3, ik]
                for it = 1:Nt
                    r1x = s.r_array[1, particle, it]
                    r1y = s.r_array[2, particle, it]
                    r1z = s.r_array[3, particle, it]
                    rk1 = kx*r1x + ky*r1y + kz*r1z
                    sinrk1, cosrk1 = sincos(rk1)
                    coskr[it, ik] = cosrk1
                    sinkr[it, ik] = sinrk1
                end
            end
        elseif s.Ndims == 2
            @turbo for ik = 1:Nk
                kx = k_array[1, ik]
                ky = k_array[2, ik]
                for it = 1:Nt
                    r1x = s.r_array[1, particle, it]
                    r1y = s.r_array[2, particle, it]
                    rk1 = kx*r1x + ky*r1y
                    sinrk1, cosrk1 = sincos(rk1)
                    coskr[it, ik] = cosrk1
                    sinkr[it, ik] = sinrk1
                end
            end
        else
            error("Only 2D and 3D simulations are supported")
        end
        for iδt in eachindex(s.dt_array)
            fkspartialidt = 0.0
            pairs_idt = s.t1_t2_pair_array[iδt]
            Npairs = size(pairs_idt, 1)
            @turbo for ipair = 1:Npairs
                t1 = pairs_idt[ipair, 1]
                t2 = pairs_idt[ipair, 2]
                for ik = 1:Nk
                    fkspartialidt += coskr[t2, ik]*coskr[t1, ik] + sinkr[t2, ik]*sinkr[t1, ik]
                end
            end
            Fks_per_particle[iδt, particle] += fkspartialidt
        end
    end

    for iδt in eachindex(s.dt_array)
        pairs_idt = s.t1_t2_pair_array[iδt]
        Npairs = size(pairs_idt, 1)
        Fks_per_particle[iδt, :] ./= Nk*Npairs
    end
    Fk = sum(Fks_per_particle; dims=2)[:] / N
    return Fk, Fks_per_particle
end



