"""
    gr3_kernel!(g_r_acc, r1, r2, task, Ntasks, N_timesteps, N1, N2, xbox_size, ybox_size, zbox_size, binsize, Nbins)

Kernel for the 3D radial distribution function.
"""
function gr3_kernel!(g_r_acc, r1, r2, task, Ntasks, N_timesteps, N1, N2, xbox_size, ybox_size, zbox_size, binsize, Nbins, verbose)
	for t = task:Ntasks:N_timesteps
        if verbose && task == 1 && t % 100 == 0
            println("Processing time step $t of $N_timesteps")
        end
		for particle1 = 1:N1
			@inbounds for particle2 = particle1+1:N2
				dx = r1[1, particle1, t] - r2[1, particle2, t]
				dy = r1[2, particle1, t] - r2[2, particle2, t]
				dz = r1[3, particle1, t] - r2[3, particle2, t]
				dx -= round(dx/xbox_size)*xbox_size
				dy -= round(dy/ybox_size)*ybox_size
				dz -= round(dz/zbox_size)*zbox_size
				@fastmath dr = sqrt(abs(dx^2+dy^2+dz^2))
				index = ceil(Int64, dr / binsize)
				index2 = ifelse(0<index<=Nbins, index, 1)
				g_r_acc[index2, task] += 2
			end
		end
	end
end

"""
    gr2_kernel!(g_r_acc, r1, r2, task, Ntasks, N_timesteps, N1, N2, xbox_size, ybox_size, binsize, Nbins)

Kernel for the 2D radial distribution function.
"""
function gr2_kernel!(g_r_acc, r1, r2, task, Ntasks, N_timesteps, N1, N2, xbox_size, ybox_size, binsize, Nbins, verbose)
	for t = task:Ntasks:N_timesteps
        if verbose && task == 1 && t % 100 == 0
            println("Processing time step $t of $N_timesteps")
        end
		for particle1 = 1:N1
			@inbounds for particle2 = particle1+1:N2
				dx = r1[1, particle1, t] - r2[1, particle2, t]
				dy = r1[2, particle1, t] - r2[2, particle2, t]
				dx -= round(dx/xbox_size)*xbox_size
				dy -= round(dy/ybox_size)*ybox_size
				@fastmath dr = sqrt(abs(dx^2+dy^2))
				index = ceil(Int64, dr / binsize)
				index2 = ifelse(0<index<=Nbins, index, 1)
				g_r_acc[index2, task] += 2
			end
		end
	end
end

"""
    find_radial_distribution_function3D(r1, r2, Nbins, box_sizes, bin_edges, Ntasks)

Calculates the 3D radial distribution function.
"""
function find_radial_distribution_function3D(r1, r2, Nbins, box_sizes, bin_edges, Ntasks, verbose)
    if verbose
        println("Calculating g(r), with $Ntasks tasks distributed over $(Threads.nthreads()) threads")
    end
    @assert length(box_sizes) == 3

    Ndim, N1, N_timesteps = size(r1)
    Ndim2, N2, N_timesteps2 = size(r2)
    @assert Ndim2 == Ndim && N_timesteps2 == N_timesteps
    rho2 = N2 / prod(box_sizes)
	#N_timesteps = div(N_timesteps,2) 

    binsize = bin_edges[2] - bin_edges[1]
    g_r_acc = zeros(Int64, Nbins, Ntasks)
    xbox_size = box_sizes[1]
    ybox_size = box_sizes[2]
    zbox_size = box_sizes[3]
    g_r_acc = zeros(Int64, Nbins, Ntasks)
        
    if Ntasks == 1
        gr3_kernel!(g_r_acc, r1, r2, 1, 1, N_timesteps, N1, N2, xbox_size, ybox_size, zbox_size, binsize, Nbins, verbose)
    else
        Threads.@threads for task = 1:Ntasks
            gr3_kernel!(g_r_acc, r1, r2, task, Ntasks, N_timesteps, N1, N2, xbox_size, ybox_size, zbox_size, binsize, Nbins, verbose)
        end
    end
	g_r_acc = sum(g_r_acc, dims=2)[:]
    g_r = float.(g_r_acc)
    for i = 1:Nbins
        volume_shell = 4π/3*(bin_edges[i + 1]^3 - bin_edges[i]^3)
        g_r[i] /= volume_shell*rho2*N_timesteps*N1
    end
    g_r[1] = 0.0
    return reshape(g_r, length(g_r))
end

"""
    find_radial_distribution_function2D(r1, r2, Nbins, box_sizes, bin_edges, Ntasks)

Calculates the 2D radial distribution function.
"""
function find_radial_distribution_function2D(r1, r2, Nbins, box_sizes, bin_edges, Ntasks, verbose)
    if verbose
        println("Calculating g(r), with $Ntasks tasks distributed over $(Threads.nthreads()) threads")
    end
    @assert length(box_sizes) == 2

    Ndim, N1, N_timesteps = size(r1)
    Ndim2, N2, N_timesteps2 = size(r2)
    @assert Ndim2 == Ndim && N_timesteps2 == N_timesteps
    rho2 = N2 / prod(box_sizes)
	#N_timesteps = div(N_timesteps,2) 

    binsize = bin_edges[2] - bin_edges[1]
    g_r_acc = zeros(Int64, Nbins, Ntasks)
    xbox_size = box_sizes[1]
    ybox_size = box_sizes[2]
    g_r_acc = zeros(Int64, Nbins, Ntasks)
    if Ntasks == 1
        gr2_kernel!(g_r_acc, r1, r2, 1, 1, N_timesteps, N1, N2, xbox_size, ybox_size, binsize, Nbins, verbose)
    else
        Threads.@threads for task = 1:Ntasks
            gr2_kernel!(g_r_acc, r1, r2, task, Ntasks, N_timesteps, N1, N2, xbox_size, ybox_size, binsize, Nbins, verbose)
        end
    end
	g_r_acc = sum(g_r_acc, dims=2)[:]
    g_r = float.(g_r_acc)
    for i = 1:Nbins
        volume_shell = π*(bin_edges[i + 1]^2 - bin_edges[i]^2)
        g_r[i] /= volume_shell*rho2*N_timesteps*N1
    end
    g_r[1] = 0.0
    return reshape(g_r, length(g_r))
end


"""
    find_radial_distribution_function(s::MultiComponentSimulation, Nbins::Int, rmax::Float64; Ntasks=1, verbose=true)

Calculates the partial radial distribution functions `g_αβ(r)` for a multi-component simulation.

The `g_αβ(r)` describes the probability of finding a particle of species `β` at a distance `r` from a particle of species `α`, relative to that of an ideal gas.

# Arguments
- `s::MultiComponentSimulation`: The simulation data.
- `Nbins::Int`: The number of bins to use for the histogram.
- `rmax::Float64`: The maximum distance to consider. If negative, it's set to half the smallest box dimension.
- `Ntasks::Int=1`: The number of parallel tasks to use for the calculation.
- `verbose::Bool=true`: If `true`, prints progress to the console.

# Returns
- `bin_centres::Vector{Float64}`: A vector of the centre points of the histogram bins.
- `g_r::Matrix{Vector{Float64}}`: A matrix of vectors, where `g_r[α, β]` is the partial `g_αβ(r)`.
"""
function find_radial_distribution_function(s::MultiComponentSimulation, Nbins, rmax; Ntasks=1, verbose=true)
    if rmax < 0.0
        rmax = minimum(s.box_sizes)/2
    end
    bin_edges = collect(LinRange(0.0, rmax, Nbins+1))
    bin_centres = [(bin_edges[i+1]+bin_edges[i])/2 for i = 1:Nbins]
    Ndims = length(s.box_sizes)
    if Ndims == 3
        return bin_centres, [find_radial_distribution_function3D(s.r_array[i], s.r_array[j], Nbins, s.box_sizes, bin_edges, Ntasks, verbose) for i=1:s.N_species, j=1:s.N_species]
    elseif Ndims == 2
        return bin_centres, [find_radial_distribution_function2D(s.r_array[i], s.r_array[j], Nbins, s.box_sizes, bin_edges, Ntasks, verbose) for i=1:s.N_species, j=1:s.N_species]
    else
        error("Not implemented yet")
    end

end

"""
    find_radial_distribution_function(s::SingleComponentSimulation, Nbins::Int, rmax::Float64; Ntasks::Int=1, verbose=true)

Calculates the radial distribution function `g(r)` for a single-component simulation.

The `g(r)` describes the probability of finding a particle at a distance `r` from another particle, relative to that of an ideal gas.

# Arguments
- `s::SingleComponentSimulation`: The simulation data.
- `Nbins::Int`: The number of bins to use for the histogram.
- `rmax::Float64`: The maximum distance to consider. If negative, it's set to half the smallest box dimension.
- `Ntasks::Int=1`: The number of parallel tasks to use for the calculation.
- `verbose::Bool=true`: If `true`, prints progress to the console.

# Returns
- `bin_centres::Vector{Float64}`: A vector of the centre points of the histogram bins.
- `g_r::Vector{Float64}`: A vector containing the values of `g(r)`.
"""
function find_radial_distribution_function(s::SingleComponentSimulation, Nbins::Int, rmax; Ntasks=1, verbose=true)
    if rmax < 0.0
        rmax = minimum(box_sizes)/2
    end
    bin_edges = collect(LinRange(0.0, rmax, Nbins+1))
    bin_centres = [(bin_edges[i+1]+bin_edges[i])/2 for i = 1:Nbins]
    Ndims = length(s.box_sizes)
    if Ndims == 3
        return bin_centres, find_radial_distribution_function3D(s.r_array, s.r_array, Nbins, s.box_sizes, bin_edges, Ntasks, verbose)
    elseif Ndims == 2
        return bin_centres, find_radial_distribution_function2D(s.r_array, s.r_array, Nbins, s.box_sizes, bin_edges, Ntasks, verbose)
    else
        error("Not implemented yet")
    end
end

"""
    find_radial_distribution_function(s::SingleComponentSimulation, dr::Float64, rmax::Float64; Ntasks::Int=1, verbose=true)

Calculates the radial distribution function `g(r)` for a single-component simulation.

This method defines the histogram bins by a fixed width `dr`.

# Arguments
- `s::SingleComponentSimulation`: The simulation data.
- `dr::Float64`: The width of the histogram bins.
- `rmax::Float64`: The maximum distance to consider. If negative, it's set to half the smallest box dimension.
- `Ntasks::Int=1`: The number of parallel tasks to use for the calculation.
- `verbose::Bool=true`: If `true`, prints progress to the console.

# Returns
- `bin_centres::Vector{Float64}`: A vector of the centre points of the histogram bins.
- `g_r::Vector{Float64}`: A vector containing the values of `g(r)`.
"""
function find_radial_distribution_function(s::SingleComponentSimulation, dr::Float64, rmax; Ntasks=1, verbose=true)
    if rmax < 0.0
        rmax = minimum(box_sizes)/2
    end
    Nbins = floor(Int, rmax/dr)
    bin_edges = range(0, dr*Nbins, length=Nbins+1)
    bin_centres = [(bin_edges[i+1]+bin_edges[i])/2 for i = 1:Nbins]
    Ndims = length(s.box_sizes)
    if Ndims == 3
        return bin_centres, find_radial_distribution_function3D(s.r_array, s.r_array, Nbins, s.box_sizes, bin_edges, Ntasks, verbose)
    elseif Ndims == 2
        return bin_centres, find_radial_distribution_function2D(s.r_array, s.r_array, Nbins, s.box_sizes, bin_edges, Ntasks, verbose)
    else
        error("Not implemented yet")
    end
end