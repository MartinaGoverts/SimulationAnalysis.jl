# TO DO: partial velocity correlations (with different forces)
# TO DO: 3D calculation of active force (this might also depend on simulation details so leave for now)
# added keyword to current mode calculation (specify which force to use)

struct SingleComponentCurrentModes <: AbstractDensityModes
    Re::Array{Float64, 2}
    Im::Array{Float64, 2}
end

# current mode: defined according to Reichman 2005 (Debets 2023 doesn't include the 1/k prefactor)
"""
    find_current_modes(s::SingleComponentSimulation, p, kspace::Kspace; verbose = true)

A function that calculates the density current modes, which are defined as:
j(k,t) = μ / |k| ∑_j (k dot Ftot) exp(i k dot r_j(t)),
where r_j is the position of particle j, μ is the particle's mobility, Ftot is the total (interaction + active) force and
k is the wavevector at which the system is probed.

# Arguments
- `s::Union{SingleComponentSimulation, SelfPropelledVoronoiSimulation}`: a single-component SelfPropelledVoronoi simulation
- `kspace::KSpace`: [...]
- `verbose=true`: [...]
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

function show(io::IO,  ::MIME"text/plain", ρkt::SingleComponentCurrentModes)
    println(io, "SingleComponentCurrentModes with real and imaginary parts of size $(size(ρkt.Re)).")
end

# calculate the total force on each particle; for all particles & time steps in the simulation
# f: interaction force (sim.F_array)
# u: particle orientations (sim.u_array)
# v0: active force strength (float)
# μ: mobility (1/friction constant)
# NOTE: only 2D so far
# Ftot = Fint + v0 / μ * orientation_vector
function calculate_total_force(f, u, v0, μ)
    @assert size(u) == size(f[1,:,:])
    Fact = similar(f)
    Fact[1,:,:] = @. v0 / μ * cos.(u)
    Fact[2,:,:] = @. v0 / μ * sin.(u)
    return Fact .+ f
end

# r: sim.r_array (positions)
# f: sim.F_array (interaction forces)
# u: sim.u_array (orientation angles)
# v0: active force strengths / species (vector)
# idea: keyword (i, a, t) - interaction force, active force or total force?
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

#####################################################################

struct MultiComponentCurrentModes <: AbstractDensityModes
    Re::Vector{Array{Float64, 2}}
    Im::Vector{Array{Float64, 2}}
end

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

function show(io::IO,  ::MIME"text/plain", ρkt::MultiComponentCurrentModes)
    println(io, "MultiComponentCurrentModes with real and imaginary parts of size $(size(ρkt.Re[1])) for $(length(ρkt.Re)) species.")
end

##########################################################################

# calculate w0 / w(infty); averaged over all particles and all time steps
function find_force_correlation(s::SelfPropelledVoronoiSimulation)
    Fint = s.F_array
    orient = s.u_array
    N = s.N; Nt = s.Nt; Ndims = s.Ndims
    v0 = s.v0; μ = s.mobility
    forces = calculate_total_force(Fint, orient, v0, μ)
    force_sum = [(μ*forces[i,:,:]).^2 for i=1:Ndims]
    return sum(sum(force_sum)) / (Ndims*Nt*N)
end

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