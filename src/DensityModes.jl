
abstract type AbstractDensityModes end

"""
    SingleComponentDensityModes

A struct to hold the density modes `ρ(k,t)` for a single-component simulation.

The density modes are defined as `ρ(k,t) = Σ_j exp(i * k ⋅ r_j(t))`, where the sum is over all particles `j`.
This struct stores the real and imaginary parts of `ρ(k,t)` for various k-vectors and time steps.

# Fields
- `Re::Matrix{Float64}`: The real part of the density modes, `Σ_j cos(k ⋅ rj)`. Dimensions are `(N_timesteps, Nk)`.
- `Im::Matrix{Float64}`: The imaginary part of the density modes, `Σ_j sin(k ⋅ rj)`. Dimensions are `(N_timesteps, Nk)`.
"""
struct SingleComponentDensityModes <: AbstractDensityModes
    Re::Array{Float64, 2}
    Im::Array{Float64, 2}
end

function show(io::IO,  ::MIME"text/plain", ρkt::SingleComponentDensityModes)
    println(io, "SingleComponentDensityModes with real and imaginary parts of size $(size(ρkt.Re)).")
end

"""
    MultiComponentDensityModes

A struct to hold the density modes `ρ_s(k,t)` for each species `s` in a multi-component simulation.

The density modes for each species are defined as `ρ_s(k,t) = Σ_j exp(i * k ⋅ r_j(t))`, where the sum is over all particles `j` of species `s`.
This struct stores the real and imaginary parts of `ρ_s(k,t)` for various k-vectors and time steps.

# Fields
- `Re::Vector{Matrix{Float64}}`: A vector where `Re[s]` is the real part of the density modes for species `s`. Each matrix has dimensions `(N_timesteps, Nk)`.
- `Im::Vector{Matrix{Float64}}`: A vector where `Im[s]` is the imaginary part of the density modes for species `s`. Each matrix has dimensions `(N_timesteps, Nk)`.
"""
struct MultiComponentDensityModes <: AbstractDensityModes
    Re::Vector{Array{Float64, 2}}
    Im::Vector{Array{Float64, 2}}
end

function show(io::IO,  ::MIME"text/plain", ρkt::MultiComponentDensityModes)
    println(io, "MultiComponentDensityModes with real and imaginary parts of size $(size(ρkt.Re[1])) for $(length(ρkt.Re)) species.")
end




"""
    find_density_modes(s::SingleComponentSimulation, kspace::KSpace; verbose=true)

Calculates the density modes `ρ(k,t) = Σ_j exp(i * k ⋅ r_j(t))` for a single-component simulation.

The calculation is performed for all k-vectors in `kspace` and for all time steps in the simulation.
This function is computationally intensive. The `verbose` option can be used to monitor progress.

# Arguments
- `s::SingleComponentSimulation`: The simulation data.
- `kspace::KSpace`: The k-space vectors for which to calculate the density modes.
- `verbose::Bool=true`: If `true`, prints progress and performance information to the console.

# Returns
- `SingleComponentDensityModes`: A struct containing the real and imaginary parts of the density modes.

See also: [`SingleComponentDensityModes`](@ref), [`find_density_modes(::Union{MultiComponentSimulation,MCSPVSimulation}, ::KSpace)`](@ref)
"""
function find_density_modes(s::Union{SingleComponentSimulation, SelfPropelledVoronoiSimulation}, kspace::KSpace; verbose=true)
    Ndim, N, N_timesteps = size(s.r_array)
    Nk = kspace.Nk
    Reρ = zeros(N_timesteps, Nk)
    Imρ = zeros(N_timesteps, Nk)
    if verbose
        println("Calculating density modes for $N particles at $N_timesteps time points for $Nk wave vectors")
        println("Memory usage: $(Base.format_bytes(2*Base.summarysize(Reρ)))")
        println("Based on 10 GFLOPS, this will take approximately $(round(Nk*s.N*N_timesteps*9/10^10, digits=1)) seconds.")
    end
    tstart = time()
    _find_density_modes!(Reρ, Imρ, s.r_array, kspace)
    if verbose
        tstop = time()
        println("Elapsed time: $(round(tstop-tstart,digits=3)) seconds")
        println("Achieved GFLOPS: $(round(Nk*s.N*N_timesteps*9/(tstop-tstart)/10^9, digits=3))")
    end
    return SingleComponentDensityModes(Reρ, Imρ)
end

"""
    find_density_modes(s::Union{MultiComponentSimulation,MCSPVSimulation}, kspace::KSpace; verbose=true)

Calculates the density modes `ρ_s(k,t) = Σ_j exp(i * k ⋅ r_j(t))` for each species `s` in a multi-component simulation.

The calculation is performed for all k-vectors in `kspace` and for all time steps in the simulation.
This function is computationally intensive. The `verbose` option can be used to monitor progress.

# Arguments
- `s::Union{MultiComponentSimulation,MCSPVSimulation}`: The simulation data.
- `kspace::KSpace`: The k-space vectors for which to calculate the density modes.
- `verbose::Bool=true`: If `true`, prints progress and performance information to the console.

# Returns
- `MultiComponentDensityModes`: A struct containing the real and imaginary parts of the density modes for each species.

See also: [`MultiComponentDensityModes`](@ref), [`find_density_modes(::SingleComponentSimulation, ::KSpace)`](@ref)
"""
function find_density_modes(s::Union{MultiComponentSimulation,MCSPVSimulation}, kspace::KSpace; verbose=true)
    N_species = length(s.r_array)
    Nk = kspace.Nk
    N_timesteps = size(s.r_array[1], 3)
    Reρ = [zeros(N_timesteps, Nk) for i = 1:N_species]
    Imρ = [zeros(N_timesteps, Nk) for i = 1:N_species]
    if verbose
        println("Calculating density modes for $(s.N) particles at $N_timesteps time points for $Nk wave vectors")
        println("Memory usage: $(Base.format_bytes(2*Base.summarysize(Reρ)))")
        println("Based on 10 GFLOPS, this will take approximately $(round(Nk*s.N*N_timesteps*9/10^10, digits=1)) seconds.")
    end
    tstart = time()
    for species in 1:N_species
        _find_density_modes!(Reρ[species], Imρ[species], s.r_array[species], kspace)
    end
    if verbose
        tstop = time()
        println("Elapsed time: $(round(tstop-tstart,digits=3)) seconds")
        println("Achieved GFLOPS: $(round(Nk*s.N*N_timesteps*9/(tstop-tstart)/10^9, digits=3))")
    end

    return MultiComponentDensityModes(Reρ, Imρ)
end




    #Benchmarks: Nk =29791 , N =1500, Nt = 1001, laptop 01/07/22

#    threads, turbo, loop: t,i,k:
#    23.1
#    batch percore, turbo, loop: t,i,k:
#    24.1
#    batch perthread, turbo, loop: t,i,k:
#    19.7
#    tturbo, loop: t,i,k:
#    27.7
#
#    threads, turbo, loop: t,k,i:
#    20.8
#    batch percore, turbo, loop: t,k,i:
#    25.1
#    batch perthread, turbo, loop: t,k,i:
#    19.7
#    tturbo, loop: t,k,i:
#    25.4
"""
    _find_density_modes!(Reρ, Imρ, r, kspace)

In-place calculation of density modes.

This is a helper function that computes `ρ(k,t)` and stores the result in pre-allocated arrays.

# Arguments
- `Reρ`: A pre-allocated `(N_timesteps, Nk)` matrix to store the real part of the density modes.
- `Imρ`: A pre-allocated `(N_timesteps, Nk)` matrix to store the imaginary part of the density modes.
- `r`: A `(Ndims, N, N_timesteps)` array of particle positions.
- `kspace`: A `KSpace` object containing the k-vectors.

# Returns
- `Nothing`: The `Reρ` and `Imρ` arrays are modified in-place.
"""
function _find_density_modes!(Reρ, Imρ, r, kspace)
    Ndim, N, N_timesteps = size(r)
    Nk = kspace.Nk
    k_array = kspace.k_array
    if Ndim == 3
        @batch per=thread for t = 1:N_timesteps
            @turbo for i_k = 1:Nk
                kx = k_array[1, i_k]
                ky = k_array[2, i_k]
                kz = k_array[3, i_k]
                Reρkt = 0.0
                Imρkt = 0.0
                for particle = 1:N 
                    rx = r[1, particle,  t]
                    ry = r[2, particle,  t]
                    rz = r[3, particle,  t]
                    kr = kx*rx + ky*ry + kz*rz
                    sinkr, coskr = sincos(kr)
                    Reρkt += coskr
                    Imρkt += sinkr
                end
                Reρ[t, i_k] = Reρkt
                Imρ[t, i_k] = Imρkt
            end
        end
    elseif Ndim == 2
        @batch per=thread for t = 1:N_timesteps
            @turbo for i_k = 1:Nk
                kx = k_array[1, i_k]
                ky = k_array[2, i_k]
                Reρkt = 0.0
                Imρkt = 0.0
                for particle = 1:N 
                    rx = r[1, particle, t]
                    ry = r[2, particle, t]
                    kr = kx*rx + ky*ry
                    sinkr, coskr = sincos(kr)
                    Reρkt += coskr
                    Imρkt += sinkr
                end
                Reρ[t, i_k] = Reρkt
                Imρ[t, i_k] = Imρkt
            end
        end
    else
        throw(ArgumentError("Only 2D and 3D simulations are supported"))
    end
    # return Reρ, Imρ
end

