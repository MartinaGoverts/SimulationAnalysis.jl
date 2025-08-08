abstract type InteractionPotential end

"""
    Weysser <: InteractionPotential

A struct for the Weysser interaction potential.

# Fields
- `ϵ::Float64`: The interaction strength.
- `δ::Float64`: The interaction range.
"""
struct Weysser <: InteractionPotential
    ϵ::Float64
    δ::Float64
end

"""
    Berthier <: InteractionPotential

A struct for the Berthier interaction potential.

# Fields
- `c0::Float64`: The c0 coefficient.
- `c2::Float64`: The c2 coefficient.
- `c4::Float64`: The c4 coefficient.
- `ζ::Float64`: The zeta coefficient.
- `σ_ratio::Float64`: The ratio of the diameters.
"""
struct Berthier <: InteractionPotential 
    c0::Float64
    c2::Float64
    c4::Float64
    ζ::Float64
    σ_ratio::Float64
end

"""
    WCA <: InteractionPotential

A struct for the WCA interaction potential.

# Fields
- `ϵ::Float64`: The interaction strength.
- `σ::Float64`: The interaction range.
"""
struct WCA <: InteractionPotential
    ϵ::Float64
    σ::Float64
end

"""
    KAWCA <: InteractionPotential

A struct for the Kob-Andersen WCA interaction potential.

# Fields
- `ϵ11::Float64`: The interaction strength between particles of type 1.
- `ϵ12::Float64`: The interaction strength between particles of type 1 and 2.
- `ϵ22::Float64`: The interaction strength between particles of type 2.
- `σ11::Float64`: The interaction range between particles of type 1.
- `σ12::Float64`: The interaction range between particles of type 1 and 2.
- `σ22::Float64`: The interaction range between particles of type 2.
"""
struct KAWCA <: InteractionPotential
    ϵ11::Float64
    ϵ12::Float64
    ϵ22::Float64
    σ11::Float64
    σ12::Float64
    σ22::Float64
end


"""
    find_mean_D(Di, Dj, U::Weysser)

Find the mean diameter for the Weysser potential.
"""
@inline function find_mean_D(Di, Dj, U::Weysser)
    return (Di + Dj)*0.5 
end

"""
    find_mean_D(Di, Dj, U::WCA)

Find the mean diameter for the WCA potential.
"""
@inline function find_mean_D(Di, Dj, U::WCA)
    return U.σ 
end

"""
    find_mean_D(Di, Dj, U::Berthier)

Find the mean diameter for the Berthier potential.
"""
@inline function find_mean_D(Di, Dj, U::Berthier)
    return (Di + Dj)*0.5 * (1.0 - U.ζ * abs(Di - Dj)) 
end




"""
    force(r_squared, mean_d_squared, U::WCA)

Calculates the force for the WCA potential.

-∇U = (48ϵ xi^12 + 24 xi^6) r^-2 vec(r), where we defined xi = dij/rij
this function returns the force without the multiplication with vec(r)
"""
@inline function force(r_squared, mean_d_squared, U::WCA)
    inv_r_squared = 1.0/r_squared
    xi2 = mean_d_squared*inv_r_squared
    xi4 = xi2*xi2
    xi6 = xi4*xi2
    xi12 = xi6*xi6
    return U.ϵ*(48*xi12 - 24*xi6) * inv_r_squared
end

"""
    force(r_squared, mean_d_squared, U::Weysser)

Calculates the force for the Weysser potential.

-∇U = 36ϵ (xi)^36 r^-2 vec(r), where we defined xi = dij/rij
this function returns the force without the multiplication with vec(r)
"""
@inline function force(r_squared, mean_d_squared, U::Weysser)
    inv_r_squared = 1.0/r_squared
    xi2 = mean_d_squared*inv_r_squared
    xi4 = xi2*xi2
    xi8 = xi4*xi4
    xi16 = xi8*xi8
    xi36 = xi16*xi16*xi4
    return 36*U.ϵ*xi36*inv_r_squared
end

"""
    force(r_squared, mean_d_squared, U::Berthier)

Calculates the force for the Berthier potential.

Soft repulsive pair potential force (before multiplication with pair vector)
"""
@inline function force(r_squared, mean_d_squared, U::Berthier)
    inv_mean_d2 = 1.0 / mean_d_squared
    invxi2 = r_squared*inv_mean_d2
    invxi4 = invxi2*invxi2
    xi14 = 1.0/(invxi4*invxi4*invxi4*invxi2)
    return -2.0 * inv_mean_d2 * (U.c2 + 2.0*U.c4*invxi2 - 6.0*xi14)
end

"""
    calculate_forces!(s, U; friction=false, cutoff=1.25)

Recalculates the total force on all particles according to the langevin equation F = -∇U - γv. This function also updates the
total potential energy in the output datastructure.
"""
function calculate_forces!(s, U::InteractionPotential; friction=false, cutoff=1.25)
    r_array = s.r_array
    F_array = s.F_array
    F_array .= 0

    _, N, N_timesteps = size(r_array)

    if typeof(s) == ContinuouslyPolydisperseSimulation
        D_array = s.D_array
    else
        D_array = ones(N)
    end

    F_array .= 0.0
    r²_cutoff = cutoff^2
    box_size = s.box_size
    @time @batch for t = 1:N_timesteps
        @inbounds for particle_i = 1:N
            xi = r_array[1, particle_i, t]
            yi = r_array[2, particle_i, t]
            zi = r_array[3, particle_i, t]
            Di = D_array[particle_i]
            fx = 0.0
            fy = 0.0
            fz = 0.0
            for particle_j = 1:particle_i-1
                xj = r_array[1, particle_j, t]
                yj = r_array[2, particle_j, t]
                zj = r_array[3, particle_j, t]
                dx = xi-xj
                dy = yi-yj
                dz = zi-zj
                dx -= round(dx/box_size)*box_size
                dy -= round(dy/box_size)*box_size
                dz -= round(dz/box_size)*box_size
                rij2 = dx^2+dy^2+dz^2
                Dj = D_array[particle_j]
                mean_d = find_mean_D(Di, Dj, U)
                mean_d_squared = mean_d*mean_d
                F = ifelse(r²_cutoff * mean_d_squared < rij2, 0.0, force(rij2, mean_d_squared, U))
                fx += F*dx
                fy += F*dy
                fz += F*dz
                F_array[1, particle_j, t] -= F*dx
                F_array[2, particle_j, t] -= F*dy
                F_array[3, particle_j, t] -= F*dz
            end
            F_array[1, particle_i, t] += fx
            F_array[2, particle_i, t] += fy
            F_array[3, particle_i, t] += fz
        end

        # Damping
        # if friction
        #     γ = 10.0
        #     @inbounds for particle_i = 1:N
        #         F_array[1, particle_i, t] += -γ*v_array[1, particle_i, t]
        #         F_array[2, particle_i, t] += -γ*v_array[2, particle_i, t]
        #         F_array[3, particle_i, t] += -γ*v_array[3, particle_i, t]
        #     end
        # end 
    end
end

"""
    force(r_squared, ϵ, σ, ::KAWCA)

Calculates the force for the KAWCA potential.

-∇U = (48ϵ xi^12 + 24 xi^6) r^-2 vec(r), where we defined xi = dij/rij
this function returns the force without the multiplication with vec(r)
"""
@inline function force(r_squared, ϵ, σ, ::KAWCA)
    inv_r_squared = 1/ r_squared
    xi2 = σ*σ*inv_r_squared
    xi4 = xi2*xi2
    xi6 = xi4*xi2
    xi12 = xi6*xi6
    return ϵ*(48*xi12 - 24*xi6) * inv_r_squared
end


"""
    get_epsilon_sigma(type1, type2, U::KAWCA)

Get the epsilon and sigma for the KAWCA potential.
"""
@inline function get_epsilon_sigma(type1, type2, U::KAWCA)
    if type1 == type2
        ϵ = ifelse(type1 == 1, U.ϵ11, U.ϵ22)
        σ = ifelse(type1 == 1, U.σ11, U.σ22)
    else 
        ϵ = U.ϵ12
        σ = U.σ12
    end
    return (ϵ, σ)
end

"""
    calculate_forces!(s, U::KAWCA; cutoff=2.0^(1.0/6.0))

Recalculates the total force on all particles in a multi-component simulation `s` according to the Kob-Andersen Weeks-Chandler-Andersen system U, and updates the force array in `s` with the calculated values.

# Arguments
- `s`: a `MultiComponentSimulation` object that contains the position and force arrays for each particle.
- `U::KAWCA`: the interaction type.
- `cutoff`: the cutoff distance for the potential energy function. Default is `2.0^(1.0/6.0)`.
"""
function calculate_forces!(s, U::KAWCA; cutoff=2.0^(1.0/6.0))
    r_array = s.r_array
    F_array = s.F_array
    N_timesteps = length(s.t_array)

    @assert typeof(s) == MultiComponentSimulation
    @assert s.box_sizes[1] == s.box_sizes[2] == s.box_sizes[3] 
    @assert s.N_species == 2
    r²_cutoff = cutoff^2
    box_size = s.box_sizes[1]
    Ns = s.N_species
    @time @batch for t = 1:N_timesteps
        @inbounds for typei = 1:Ns
            r_array_i = r_array[typei]
              for particle_i = 1:s.N_particles_per_species[typei]
                xi = r_array_i[1, particle_i, t]
                yi = r_array_i[2, particle_i, t]
                zi = r_array_i[3, particle_i, t]
                fx = 0.0
                fy = 0.0
                fz = 0.0
                for typej = 1:Ns
                    r_array_j = r_array[typej]

                    for particle_j = 1:s.N_particles_per_species[typej]
                        if typei == typej && particle_i == particle_j
                            continue
                        end
                        xj = r_array_j[1, particle_j, t]
                        yj = r_array_j[2, particle_j, t]
                        zj = r_array_j[3, particle_j, t]
                        dx = xi-xj
                        dy = yi-yj
                        dz = zi-zj
                        dx -= round(dx/box_size)*box_size
                        dy -= round(dy/box_size)*box_size
                        dz -= round(dz/box_size)*box_size
                        rij2 = dx^2+dy^2+dz^2
                        (ϵ, σ) = get_epsilon_sigma(typei, typej, U)

                        F = ifelse(r²_cutoff * σ^2 < rij2, 0.0, force(rij2, ϵ, σ, U))

                        fx += F*dx
                        fy += F*dy
                        fz += F*dz
                    end
                end
                F_array[typei][1, particle_i, t] = fx
                F_array[typei][2, particle_i, t] = fy
                F_array[typei][3, particle_i, t] = fz
            end
        end
    end
end