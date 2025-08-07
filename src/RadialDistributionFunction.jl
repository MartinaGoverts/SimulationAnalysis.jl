"""
    gr3_kernel!(g_r_acc, r1, r2, task, Ntasks, N_timesteps, N1, N2, xbox_size, ybox_size, zbox_size, binsize, Nbins)

Kernel for the 3D radial distribution function.
"""
function gr3_kernel!(g_r_acc, r1, r2, task, Ntasks, N_timesteps, N1, N2, xbox_size, ybox_size, zbox_size, binsize, Nbins)
	for t = task:Ntasks:N_timesteps
		if task == 1
			@show t, N_timesteps
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
function gr2_kernel!(g_r_acc, r1, r2, task, Ntasks, N_timesteps, N1, N2, xbox_size, ybox_size, binsize, Nbins)
	for t = task:Ntasks:N_timesteps
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
function find_radial_distribution_function3D(r1, r2, Nbins, box_sizes, bin_edges, Ntasks)
    println("Calculating g(r), with $Ntasks tasks distributed over $(Threads.nthreads()) threads")
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
    @time Threads.@threads for task = 1:Ntasks
		gr3_kernel!(g_r_acc, r1, r2, task, Ntasks, N_timesteps, N1, N2, xbox_size, ybox_size, zbox_size, binsize, Nbins)
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
function find_radial_distribution_function2D(r1, r2, Nbins, box_sizes, bin_edges, Ntasks)
    println("Calculating g(r), with $Ntasks tasks distributed over $(Ntasks) threads")
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
    Threads.@threads for task = 1:Ntasks
		gr2_kernel!(g_r_acc, r1, r2, task, Ntasks, N_timesteps, N1, N2, xbox_size, ybox_size, binsize, Nbins)
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
    find_radial_distribution_function2D2(r1, r2, Nbins, box_sizes, bin_edges)

Calculates the 2D radial distribution function.
"""
function find_radial_distribution_function2D2(r1, r2, Nbins, box_sizes, bin_edges)
    println("Calculating g(r)")
    @assert length(box_sizes) == 2
    Ndim, N1, N_timesteps = size(r1)
    Ndim2, N2, N_timesteps2 = size(r2)
    @assert Ndim2 == Ndim && N_timesteps2 == N_timesteps
    rho2 = N2 / prod(box_sizes)

    binsize = bin_edges[2] - bin_edges[1]
    g_r_acc = zeros(Int64, Nbins)
    xbox_size = box_sizes[1]
    ybox_size = box_sizes[2]

    @time for t = 1:N_timesteps
        if t % 100 == 0
            @show t, N_timesteps
        end
        for particle1 = 1:N1
            for particle2 = 1:N2
                dx = r1[1, particle1, t] - r2[1, particle2, t]
                dy = r1[2, particle1, t] - r2[2, particle2, t]
                dx -= round(dx/xbox_size)*xbox_size
                dy -= round(dy/ybox_size)*ybox_size
                @fastmath dr = sqrt(abs(dx^2+dy^2))
                index = ceil(Int64, dr / binsize)
                index2 = ifelse(0 < index <= Nbins, index, 1)
                g_r_acc[index2] += 1
            end
        end
    end
    g_r = float.(g_r_acc)
    for i = 1:Nbins
        volume_shell = π*(bin_edges[i + 1]^2 - bin_edges[i]^2)
        g_r[i] /= volume_shell*rho2*N_timesteps*N1
    end
    return reshape(g_r, length(g_r))
end

"""
    find_radial_distribution_function(s::MultiComponentSimulation, Nbins, rmax)

Calculates the radial distribution function for a multi-component simulation.
"""
function find_radial_distribution_function(s::MultiComponentSimulation, Nbins, rmax)
    if rmax < 0.0
        rmax = minimum(s.box_sizes)/2
    end
    bin_edges = collect(LinRange(0.0, rmax, Nbins+1))
    bin_centres = [(bin_edges[i+1]+bin_edges[i])/2 for i = 1:Nbins]
    Ndims = length(s.box_sizes)
    if Ndims == 3
        return bin_centres, [find_radial_distribution_function3D(s.r_array[i], s.r_array[j], Nbins, s.box_sizes, bin_edges) for i=1:s.N_species, j=1:s.N_species]
    elseif Ndims == 2
        return bin_centres, [find_radial_distribution_function2D(s.r_array[i], s.r_array[j], Nbins, s.box_sizes, bin_edges) for i=1:s.N_species, j=1:s.N_species]
    else
        error("Not implemented yet")
    end

end

"""
    find_radial_distribution_function(s::SingleComponentSimulation, Nbins::Int, rmax; Ntasks=1)

Calculates the radial distribution function for a single-component simulation.
"""
function find_radial_distribution_function(s::SingleComponentSimulation, Nbins::Int, rmax; Ntasks=1)
    if rmax < 0.0
        rmax = minimum(box_sizes)/2
    end
    bin_edges = collect(LinRange(0.0, rmax, Nbins+1))
    bin_centres = [(bin_edges[i+1]+bin_edges[i])/2 for i = 1:Nbins]
    Ndims = length(s.box_sizes)
    if Ndims == 3
        return bin_centres, find_radial_distribution_function3D(s.r_array, s.r_array, Nbins, s.box_sizes, bin_edges, Ntasks)
    elseif Ndims == 2
        return bin_centres, find_radial_distribution_function2D(s.r_array, s.r_array, Nbins, s.box_sizes, bin_edges, Ntasks)
    else
        error("Not implemented yet")
    end
end

"""
    find_radial_distribution_function(s::SingleComponentSimulation, dr::Float64, rmax; Ntasks=1)

Calculates the radial distribution function for a single-component simulation.
"""
function find_radial_distribution_function(s::SingleComponentSimulation, dr::Float64, rmax; Ntasks=1)
    if rmax < 0.0
        rmax = minimum(box_sizes)/2
    end
    Nbins = floor(Int, rmax/dr)
    bin_edges = range(0, dr*Nbins, length=Nbins+1)
    bin_centres = [(bin_edges[i+1]+bin_edges[i])/2 for i = 1:Nbins]
    Ndims = length(s.box_sizes)
    if Ndims == 3
        return bin_centres, find_radial_distribution_function3D(s.r_array, s.r_array, Nbins, s.box_sizes, bin_edges, Ntasks)
    elseif Ndims == 2
        return bin_centres, find_radial_distribution_function2D(s.r_array, s.r_array, Nbins, s.box_sizes, bin_edges, Ntasks)
    else
        error("Not implemented yet")
    end
end