
abstract type AbstractDensityModes end

"""
    SingleComponentDensityModes

A struct to hold the density modes of a single-component simulation.

# Fields
- `Re::Array{Float64, 2}`: The real part of the density modes.
- `Im::Array{Float64, 2}`: The imaginary part of the density modes.
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

A struct to hold the density modes of a multi-component simulation.

# Fields
- `Re::Vector{Array{Float64, 2}}`: The real part of the density modes for each species.
- `Im::Vector{Array{Float64, 2}}`: The imaginary part of the density modes for each species.
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

Calculates the density modes for a single-component simulation.

# Arguments
- `s::SingleComponentSimulation`: The simulation.
- `kspace::KSpace`: The k-space.
- `verbose::Bool=true`: Whether to print verbose output.

# Returns
- `SingleComponentDensityModes`: The density modes.
"""
function find_density_modes(s::SingleComponentSimulation, kspace::KSpace; verbose=true)
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
    find_density_modes(s::MultiComponentSimulation, kspace::KSpace; verbose=true)

Calculates the density modes for a multi-component simulation.

# Arguments
- `s::MultiComponentSimulation`: The simulation.
- `kspace::KSpace`: The k-space.
- `verbose::Bool=true`: Whether to print verbose output.

# Returns
- `MultiComponentDensityModes`: The density modes.
"""
function find_density_modes(s::MultiComponentSimulation, kspace::KSpace; verbose=true)
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

Calculate the density modes for a given simulation `r` and wave vector space `kspace`.

Reρ and Imρ arrays are preallocated arrays with size (N_timesteps, Nk) for each species of the simulation.

# Arguments
- `Reρ`: preallocated array to store the real part of the density modes.
- `Imρ`: preallocated array to store the imaginary part of the density modes.
- `r`: array containing the positions of the particles at each time step. It has shape (Ndim, N, N_timesteps), where
             Ndim is the number of dimensions, N is the number of particles and N_timesteps is the number of time steps.
- `kspace`: the wave vector space used to calculate the density modes.

# Returns
- `Nothing`: The `Reρ` and `Imρ` arrays are updated in-place.
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

