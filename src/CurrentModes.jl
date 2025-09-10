# TO DO: different forces --> include type of force in output / print text?
# TO DO: 3D calculation of active force (this might also depend on simulation details so leave for now)

"""
    SingleComponentCurrentModes

A struct which holds the current modes of a single-component simulation.
    
The current modes are defined as: j(k,t) = μ / |k| ∑_j (k ⋅ F_j(t)) exp(i k ⋅ r_j(t)), where μ is the
mobility, k is the wavevector, r_j(t) is the position of particle j at time t, and F_j(t) is the force on 
particle j.

# Fields
- `Re::Matrix{Float64}`: Real part of the current modes, with dimensions  `(N_timesteps, Nk)`.
- `Im::Matrix{Float64}`: Imaginary part of the current modes, `Σ_j sin(k ⋅ rj)`, with dimensions `(N_timesteps, Nk)`.
"""
struct SingleComponentCurrentModes <: AbstractDensityModes
    Re::Array{Float64, 2}
    Im::Array{Float64, 2}
end

"""
    MultiComponentCurrentModes

A struct which holds the current modes of a multi-component simulation.
    
The current modes are defined as: j{α}(k,t) = μ / |k| ∑_j (k ⋅ F{α}_j(t)) exp(i k ⋅ r{α}_j(t)),
where {α} denotes the subspecies, and the sum over j now runs over the number of particles of subspecies α.

# Fields
- `Re::Matrix{Float64}`: Real part of the current modes, with dimensions  `(N_timesteps, Nk)`.
- `Im::Matrix{Float64}`: Imaginary part of the current modes, `Σ_j sin(k ⋅ rj)`, with dimensions `(N_timesteps, Nk)`.
"""
struct MultiComponentCurrentModes <: AbstractDensityModes
    Re::Vector{Array{Float64, 2}}
    Im::Vector{Array{Float64, 2}}
end


# if we include partial velocity correlations, update the docstring here
"""
    find_current_modes(s::Union{SingleComponentSimulation, SelfPropelledVoronoiSimulation}, kspace::KSpace; verbose=true)

A function which calculates the current modes, which are defined as: j(k,t) = μ / |k| ∑_j (k ⋅ F_j(t)) exp(i k ⋅ r_j(t)).
Here r_j(t) is the position of particle j, μ is the particle's mobility, F_j(t) is the force on particle j and
k is the wavevector at which the system is probed.

The total force (interaction + active) on every particle is calculated. The current modes are stored for every timestep in simulation `s`, and
for every wavevector in `kspace`.

# Arguments
- `s::Union{SingleComponentSimulation, SelfPropelledVoronoiSimulation}`: A single-component SelfPropelledVoronoi simulation
- `kspace::KSpace`: A Kspace object, which contains the wavevectors at which to evaluate the current modes.
- `verbose=true`: If `true`, prints a performance overview.

# Returns
A SingleComponentCurrentModes object, with fields `Re` and `Im` that contain the real and imaginary parts of the current modes
at all times and for all wavevectors.
"""
function find_current_modes(s::Union{SingleComponentSimulation, SelfPropelledVoronoiSimulation}, kspace::KSpace; verbose=true)
    Ndim, N, N_timesteps = size(s.r_array)
    Nk = kspace.Nk

    Rej = zeros(N_timesteps, Nk)
    Imj = zeros(N_timesteps, Nk)

    if verbose
        println("Calculating current modes for $N particles at $N_timesteps time points for $Nk wave vectors")
        println("Memory usage: $(Base.format_bytes(2*Base.summarysize(Rej)))")
        println("Based on 10 GFLOPS, this will take approximately $(round(Nk*s.N*N_timesteps*9/10^10, digits=1)) seconds.")
    end
    tstart = time()

    _find_current_modes!(Rej, Imj, s.r_array, s.F_array, s.u_array, s.v0, s.mobility, kspace)

    if verbose
        tstop = time()
        println("Elapsed time: $(round(tstop-tstart,digits=3)) seconds")
        println("Achieved GFLOPS: $(round(Nk*s.N*N_timesteps*9/(tstop-tstart)/10^9, digits=3))")
    end

    return SingleComponentCurrentModes(Rej, Imj)
end


# NOTE: this only works in 2D so far. (Change w(k) so it only accepts 2D simulations?)
"""
    calculate_total_force(f, u, v0, μ)

Calculates the total force on every particle for all timesteps.
    
The total force is defined as: Ftot_j = Fint_j + v0 / μ n_j,
where Fint_j is the interaction force, and n_j = (cosθ_j, sinθ_j) is the
orientation vector of the active force. Both forces are time-dependent.

# Arguments
- `f::Array{Float64, 3}`: Contains the interaction forces of every particle. Dimensions are (2, N_particles, N_timesteps) 
- `u::Matrix{Float64}`: Contains the orientation vectors θ_j of every particle. Dimensions are (N_particles, N_timesteps)
- `v0::Float64`: Strength of the active force.
- `μ::Float64`: The mobility (1 / friction constant) of the particles.
"""
function calculate_total_force(f::Array{Float64, 3}, u::Matrix{Float64}, v0::Float64, μ::Float64)
    @assert size(u) == size(f[1,:,:])
    Fact = similar(f)
    Fact[1,:,:] = @. v0 / μ * cos.(u)
    Fact[2,:,:] = @. v0 / μ * sin.(u)
    return Fact .+ f
end

# NOTE: calculate_total_force() only 2D! (Edit later)
"""
    _find_current_modes!(Rej, Imj, r, f, u, v0, μ, kspace; force_keyword="t")

Calculates the current modes by modifying `Rej` and `Imj` in-place. Note that currently only
two-dimensional force calculations are implemented!

# Arguments
- `Rej`: The real part of the CurrentModes object, with dimensions (N_timesteps, Nk)
- `Imj`: The imaginary part of the CurrentModes object, with dimensions (N_timesteps, Nk)
- `r`: An array containing the particle positions, with dimensions (Ndims, N, N_timesteps)
- `f`: An array containing the interaction forces, with dimensions (Ndims, N, N_timesteps)
- `u`: A matrix containing the active forces, with dimensions (N, N_timesteps)
- `v0`: Active force strength.
- `μ`: Mobility.
- `kspace`: An object containing the wavevectors at which to evaluate the current modes.
- `force_keyword=t`: Specifies whether to calculate using the total force (`t`), the interaction force (`i`) or the 
        active force (`a`). The default is `t`. When a different keyword is given, the total force is calculated.
"""
function _find_current_modes!(Rej, Imj, r, f, u, v0, μ, kspace; force_keyword="t")
    Ndim, N, N_timesteps = size(r)
    Nk = kspace.Nk
    k_array = kspace.k_array
    klengths = kspace.k_lengths
    Ftot = calculate_total_force(f, u, v0, μ);

    if force_keyword == "t" # use total force
        F = Ftot
    elseif force_keyword == "i" # use interaction force
        F = f
    elseif force_keyword == "a"  # active force
        F = Ftot .- f
    else
        @warn "Received unsupported force keyword. Using total force for current modes."
        F = Ftot
    end

    if Ndim == 3
        @batch per=thread for t = 1:N_timesteps
            @turbo for i_k = 1:Nk
                kx = k_array[1, i_k]
                ky = k_array[2, i_k]
                kz = k_array[3, i_k]
                kmag = klengths[i_k]

                Rejkt = 0.0
                Imjkt = 0.0
                
                for particle = 1:N 
                    rx = r[1, particle, t]
                    ry = r[2, particle, t]
                    rz = r[3, particle, t]
                    fx = F[1, particle, t]
                    fy = F[2, particle, t]
                    fz = F[3, particle, t]

                    kr = kx*rx + ky*ry + kz*rz
                    kf = kx*fx + ky*fy + kz*fz

                    sinkr, coskr = sincos(kr)
                    Rejkt += μ * kf * coskr / kmag
                    Imjkt += μ * kf * sinkr / kmag
                end
                Rej[t, i_k] = Rejkt
                Imj[t, i_k] = Imjkt
            end
        end
    elseif Ndim == 2
        @batch per=thread for t = 1:N_timesteps
            @turbo for i_k = 1:Nk
                kx = k_array[1, i_k]
                ky = k_array[2, i_k]
                kmag = klengths[i_k]

                Rejkt = 0.0
                Imjkt = 0.0

                for particle = 1:N 
                    rx = r[1, particle, t]
                    ry = r[2, particle, t]
                    fx = F[1, particle, t]
                    fy = F[2, particle, t]

                    kr = kx*rx + ky*ry
                    kf = kx*fx + ky*fy

                    sinkr, coskr = sincos(kr)
                    Rejkt += kf * coskr / kmag
                    Imjkt += kf * sinkr / kmag
                end
                Rej[t, i_k] = Rejkt
                Imj[t, i_k] = Imjkt
            end
        end
    else
        throw(ArgumentError("Only 2D and 3D simulations are supported"))
    end
end


"""
    find_current_modes(s::Union{MultiComponentSimulation,MCSPVSimulation}, kspace::KSpace; verbose=true)

Calculates the current modes for all species in a multi-component simulation. The current modes are defined as:
j{α}(k,t) = μ / |k| ∑_j (k ⋅ F{α}_j(t)) exp(i k ⋅ r{α}_j(t)),
where {α} denotes the particle species for which to evaluate the current modes. 

The total force (interaction + active) on every particle is calculated. The current modes are stored for every timestep in simulation `s`, and
for every wavevector in `kspace`.

# Arguments
- `s::Union{MultiComponentSimulation, MCSPVSimulation}`: A multi-component SPV simulation object.
- `kspace::KSpace`: A `KSpace` object containing the wavevectors at which to evaluate the current modes.
- `verbose=true`: If `true`, prints a performance overview.

# Returns
A SingleComponentCurrentModes object, with fields `Re` and `Im` that contain the real and imaginary parts of the current modes
at all times and for all wavevectors.
"""
function find_current_modes(s::Union{MultiComponentSimulation,MCSPVSimulation}, kspace::KSpace; verbose=true)
    N_species = length(s.r_array)
    Nk = kspace.Nk
    N_timesteps = size(s.r_array[1], 3)
    klengths = kspace.k_lengths

    Rej = [zeros(N_timesteps, Nk) for i = 1:N_species]
    Imj = [zeros(N_timesteps, Nk) for i = 1:N_species]
    if verbose
        println("Calculating density modes for $(s.N) particles at $N_timesteps time points for $Nk wave vectors")
        println("Memory usage: $(Base.format_bytes(2*Base.summarysize(Rej)))")
        println("Based on 10 GFLOPS, this will take approximately $(round(Nk*s.N*N_timesteps*9/10^10, digits=1)) seconds.")
    end
    tstart = time()

    for species in 1:N_species
        _find_current_modes!(Rej[species], Imj[species], s.r_array[species], s.F_array[species], s.u_array[species], s.v0[species], s.mobility[species], kspace)
    end

    if verbose
        tstop = time()
        println("Elapsed time: $(round(tstop-tstart,digits=3)) seconds")
        println("Achieved GFLOPS: $(round(Nk*s.N*N_timesteps*9/(tstop-tstart)/10^9, digits=3))")
    end

    return MultiComponentCurrentModes(Rej, Imj)
end


"""
   find_force_correlation(s::SelfPropelledVoronoiSimulation)
   
Calculates the average magnitude squared of the total force on a particle, scaled by the mobility. This
corresponds with local velocity correlations (for k → ∞).

The correlation is defined as: 1/(d N μ^2) ∑_j < (Ftot_j)^2 >,

where d is the dimension, N is the number of particles, μ is the mobility and Ftot_j is the total force
on particle j. The sum runs over all particles, and the average <...> is performed over all timesteps in
the simulation.

# Arguments:
- `s::SelfPropelledVoronoiSimulation`: A SPV simulation object.

# Returns:
A Float64 (The average force squared)
"""
function find_force_correlation(s::SelfPropelledVoronoiSimulation)
    Fint = s.F_array
    orient = s.u_array
    N = s.N; Nt = s.Nt; Ndims = s.Ndims
    v0 = s.v0; μ = s.mobility
    forces = calculate_total_force(Fint, orient, v0, μ)
    force_sum = [(μ*forces[i,:,:]).^2 for i=1:Ndims]
    return sum(sum(force_sum)) / (Ndims*Nt*N)
end

"""
   find_force_correlation(s::MCSPVSimulation)
   
Calculates the average magnitude squared of the total force on a particle, scaled by the mobility. This
corresponds with local velocity correlations (for k → ∞). This is done for all species in the simulation.

The correlation is defined as: δ{αβ} /(d N μ{α}^2) ∑_j < (Ftot{α}_j)^2 >,

where d is the dimension, N is the number of particles, μ is the mobility and Ftot_j is the total force
on particle j. δ{αβ} is a Kronecker delta, and {α} means that the particles of subspecies α are considered.
The average is performed over all timesteps, for all species in the simulation.

# Arguments:
- `s::MCSPV`: A multi-component SPV simulation object.

# Returns:
A Matrix of size (N_species x N_species), containing the average squared force of each species.
"""
function find_force_correlation(s::MCSPVSimulation)
    N_species = s.N_species; N = s.N; Nt = s.Nt
    Ndims = s.Ndims
    Fint = s.F_array
    orient = s.u_array
    v0 = s.v0; μ = s.mobility
    w0 = zeros(N_species, N_species)

    for spec in 1:N_species
        var = calculate_total_force(Fint[spec], orient[spec], v0[spec], μ[spec])
        var_sum = [(μ[spec]*var[i,:,:]).^2 for i=1:Ndims]
        w0[spec, spec] = sum(sum(var_sum))
    end
    return w0 ./ (Ndims * N * Nt)
end

function show(io::IO,  ::MIME"text/plain", ρkt::SingleComponentCurrentModes)
    println(io, "SingleComponentCurrentModes with real and imaginary parts of size $(size(ρkt.Re)).")
end

function show(io::IO,  ::MIME"text/plain", ρkt::MultiComponentCurrentModes)
    println(io, "MultiComponentCurrentModes with real and imaginary parts of size $(size(ρkt.Re[1])) for $(length(ρkt.Re)) species.")
end