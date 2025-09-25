
"""
    read_SPV_simulation(traj, params; dt_array=nothing, t1_t2_pair_array=nothing)

reads the SPV simulation data from the specified file. If `dt_array` and `t1_t2_pair_array` are not provided, they will be automatically determined from the data.
see SelfPropelledVoronoi.jl

# Arguments
    - `traj`: The trajectory data for SelfPropelledVoronoi.jl
    - `params`: The parameters for the simulation from SelfPropelledVoronoi.jl
    - `dt_array`: The time step array.
    - `t1_t2_pair_array`: The time pairs array.
    - `COM_correction::Bool=true`: Whether to apply a correction for the collective displacement.
       We assume that all particles have the same mass.
    - `original::Bool=false`: Whether to reconstruct the original trajectories.

# Returns
    - `SelfPropelledVoronoiSimulation`: A `SelfPropelledVoronoiSimulation` object.
"""
function read_SPV_simulation(traj, params; dt_array=nothing, t1_t2_pair_array=nothing, COM_correction::Bool=true, original::Bool=false)
    # traj, params = SelfPropelledVoronoi.load_trajectory(filenamefull)
    # traj contains fields:
    # positions_trajectory::Vector{Vector{SVector{2, Float64}}}
    # orientations_trajectory::Vector{Vector{Float64}}
    # forces_trajectory::Vector{Vector{SVector{2, Float64}}}
    # potential_energy_trajectory::Vector{Float64}
    # areas_trajectory::Vector{Vector{Float64}}
    # perimeters_trajectory::Vector{Vector{Float64}}
    # steps_saved::Vector{Int64}
    r = traj.positions_trajectory |> stack |> stack
    u = traj.orientations_trajectory |> stack 
    f = traj.forces_trajectory |> stack |> stack
    Epot = traj.potential_energy_trajectory 
    A = traj.areas_trajectory |> stack
    P = traj.perimeters_trajectory |> stack
    steps_saved = traj.steps_saved
    v0 = params.particles.active_force_strengths[1]
    mobility = 1.0 / params.frictionconstant
    @assert all([params.particles.active_force_strengths[i] == v0 for i=1:params.N]) "Particles belonging to one species should have the same active force strength!"

    dt = params.dt
    dims = size(r, 1)
    t = steps_saved * dt
    box_sizes = params.box.box_sizes
    if dt_array === nothing && t1_t2_pair_array === nothing
        dt_array, t1_t2_pair_array = find_allowed_dt_array(steps_saved::Vector{Int})
        dt_array = dt_array * dt
    end
    N = size(r, 2)

    if COM_correction
        r = COM_correction_function(r,box_sizes, original)
    elseif original
        r = find_original_trajectories(r, box_sizes)
    end
   

    return SelfPropelledVoronoiSimulation(
            N, 
            dims,
            length(t),
            dt,
            v0,
            mobility,
            r,
            u,
            f,
            P,
            A,
            Epot,
            t,
            box_sizes,
            dt_array,
            t1_t2_pair_array,
            "unknown"
        )
end

"""
    read_SPV_simulation_multicomponent(traj, params, species::Vector{Int}; dt_array=nothing, t1_t2_pair_array=nothing,  COM_correction::Bool=true, original::Bool=false)

reads the SPV simulation data from the specified file. If `dt_array` and `t1_t2_pair_array` are not provided, they will be automatically determined from the data. 
Returns a multicomponent simulation object, where the species are specified by the `species` vector.
see SelfPropelledVoronoi.jl

# Arguments
    - `traj`: The trajectory data for SelfPropelledVoronoi.jl
    - `params`: The parameters for the simulation from SelfPropelledVoronoi.jl
    - `species::Vector{Int}`: A vector specifying the species of the particles.
    - `dt_array`: The time step array.
    - `t1_t2_pair_array`: The time pairs array.
    - `COM_correction::Bool=true`: Whether to apply a correction for the collective displacement.
       We assume that all particles have the same mass.
    - `original::Bool=false`: Whether to reconstruct the original trajectories.


# Returns
    - `SelfPropelledVoronoiSimulation`: A `MulticomponentSelfPropelledVoronoiSimulation` object.
"""
function read_SPV_simulation_multicomponent(traj, params, species::Vector{Int}; dt_array=nothing, t1_t2_pair_array=nothing, COM_correction::Bool=true, original::Bool=false)
    # traj, params = SelfPropelledVoronoi.load_trajectory(filenamefull)
    # traj contains fields:
    # positions_trajectory::Vector{Vector{SVector{2, Float64}}}
    # orientations_trajectory::Vector{Vector{Float64}}
    # forces_trajectory::Vector{Vector{SVector{2, Float64}}}
    # potential_energy_trajectory::Vector{Float64}
    # areas_trajectory::Vector{Vector{Float64}}
    # perimeters_trajectory::Vector{Vector{Float64}}
    # steps_saved::Vector{Int64}
    r = traj.positions_trajectory |> stack |> stack
    u = traj.orientations_trajectory |> stack 
    f = traj.forces_trajectory |> stack |> stack
    Epot = traj.potential_energy_trajectory 
    A = traj.areas_trajectory |> stack
    P = traj.perimeters_trajectory |> stack
    steps_saved = traj.steps_saved
    v0 = params.particles.active_force_strengths

    dt = params.dt
    dims = size(r, 1)  # dims = length(r[1][1])
    t = steps_saved * dt
    box_sizes = params.box.box_sizes
    if dt_array === nothing && t1_t2_pair_array === nothing
        dt_array, t1_t2_pair_array = find_allowed_dt_array(steps_saved::Vector{Int})
        dt_array = dt_array * dt
    end
    N = size(r, 2)

    if COM_correction
        r = COM_correction_function(r, box_sizes, original)
    elseif original
        r = find_original_trajectories(r, box_sizes)
    end


    # group by species

    # ensure species is a vector of integers, containing consequtive integers starting from 1
    if minimum(species) != 1 || maximum(species) != length(unique(species))
        error("Species vector must contain consecutive integers starting from 1.")
    end
    N_species = maximum(species)
    N_particles_per_species = [sum(species .== i) for i in 1:N_species]
    rvec = Array{Array{Float64, 3}, 1}(undef, N_species)
    uvec = Array{Array{Float64, 2}, 1}(undef, N_species)
    fvec = Array{Array{Float64, 3}, 1}(undef, N_species)
    Avec = Array{Array{Float64, 2}, 1}(undef, N_species)
    Pvec = Array{Array{Float64, 2}, 1}(undef, N_species)
    vvec = Array{Float64, 1}(undef, N_species)
    μvec = Array{Float64, 1}(undef, N_species)  
    # friction constant in the SPV model is the same for all species, but here we treat it "per species" for generality
    Epotvec = Array{Float64, 1}(undef, N_species)
    for i in 1:N_species
        indices = findall(species .== i)
        rvec[i] = r[:, indices, :]
        uvec[i] = u[indices, :]
        fvec[i] = f[:, indices, :]
        Avec[i] = A[indices, :]
        Pvec[i] = P[indices, :]
        # Epotvec[i] = Epot[indices]   # this is not saved per particle / species
        @assert all( [v0[indices][i] == v0[indices][1] for i=eachindex(indices)] ) "Particles belonging to one species should have the same active force strength!"
        vvec[i] = v0[indices][1]
        μvec[i] = 1.0 / params.frictionconstant  # THIS IS THE SAME FOR EVERY SPECIES! (can change later if needed)
    end
    Epotvec = Epot

    return MCSPVSimulation(
            N, 
            dims,
            N_species,
            length(t),
            dt,
            vvec,
            μvec,
            N_particles_per_species,
            rvec,
            uvec,
            fvec,
            Pvec,
            Avec,
            Epotvec,
            t,
            box_sizes,
            dt_array,
            t1_t2_pair_array,
            "unknown"  # filenamefull
        )
end





"""
    read_WCA_simulation(filenamefull, dt; maxt=-1, every=1, original=false)

Read a simulation of WCA particles done with LAMMPS. in the data file, the columns correspond to the following:

1. Atom ID
2. x position
3. y position
4. z position
5. x force
6. y force
7. z force

# Arguments
- `filenamefull::String`: Path to the simulation file.
- `dt::Float64`: Time step.
- `maxt::Int=-1`: Maximum number of time steps to read.
- `every::Int=1`: Read every `every`-th time step.
- `original::Bool=false`: Whether to reconstruct the original trajectories.

# Returns
- `SingleComponentSimulation`: A `SingleComponentSimulation` object.
"""
function read_WCA_simulation(filenamefull, dt; maxt=-1, every=1, original=false)
    println("Reading data file")
    f = open(filenamefull)
    r = Array{Array{Float64, 2}, 1}()
    F = Array{Array{Float64, 2}, 1}()
    t = Vector{Float64}()
    iter = eachline(f)
    box_size = 0
    timestep = -1
    for line in iter
        if line == "ITEM: TIMESTEP"
            timestep += 1

            if length(r) >= maxt && maxt > 0
                break
            end
            push!(t, parse(Int64, iterate(iter)[1]))
            if timestep % every != 0
                continue
            end
            if length(r) % 100 == 0
                println(length(r), " snapshots found.")
            end
            _ = iterate(iter) # ITEM:NUMBER OF ATOMS
            N = parse(Int64, iterate(iter)[1])
            _ = iterate(iter) # ITEM: BOX BOUNDS
            box_size = parse(Float64, split(iterate(iter)[1])[2])
            _ = iterate(iter)
            _ = iterate(iter)
            _ = iterate(iter) # ITEM: ATOMS
            rnew = zeros(N, 3)
            fnew = zeros(N, 3)
            for i = 1:N
                index = Parsers.parse(Int64, readuntil(f, ' '))
                # if timestep == 0
                #     readuntil(f, ' ')
                # end

                for j = 1:3
                    rnew[index, j] = Parsers.parse(Float64, readuntil(f, ' '))
                end
                for j = 1:2
                    fnew[index, j] = Parsers.parse(Float64, readuntil(f, ' '))  
                end
                fnew[index, 3] = Parsers.parse(Float64, readline(f))  
                # readline(f)
            end
            push!(r, rnew)
            push!(F, fnew)
        end
    end
    close(f)
    box_sizes = [box_size, box_size, box_size]
    t .-= 1
    println(length(r), " timesteps found with ", size(r[1])[1], " atoms.")
    r = reshape_data(r)
    remap_positions!(r, box_sizes)
    if original
        r = find_original_trajectories(r, box_sizes)
    end
    F = reshape_data(F)
    D = ones(N)
    N = size(F, 2)
    v = zeros(size(F)...)
    dt_arr, t1_t2_pair_array = find_allowed_t1_t2_pair_array_quasilog(t;doublefactor=200)
    # dt_arr = [t[i]-t[1] for i in 1:length(t)]
    # t1_t2_pair_array = [[1;; i] for i in 1:length(t)]
    s = SingleComponentSimulation(N, 3, length(t), dt, r, v, F, D, t*dt , box_sizes, dt_arr*dt , t1_t2_pair_array, filenamefull)
    return s
end

"""
    read_Newtonian_KAWCA_simulation(filenamefull, dt; maxt=-1, every=1, original=false)

Read a Newtonian simulation of a Kob-Andersen WCA mixture. The file was generated with LAMMPS. 
The columns in the data file correspond to the following:

1. Atom ID
2. Atom type
3. x position
4. y position
5. z position
6. x velocities
7. y velocities
8. z velocities

# Arguments
- `filenamefull::String`: Path to the simulation file.
- `dt::Float64`: Time step.
- `maxt::Int=-1`: Maximum number of time steps to read.
- `every::Int=1`: Read every `every`-th time step.
- `original::Bool=false`: Whether to reconstruct the original trajectories.

# Returns
- `MultiComponentSimulation`: A `MultiComponentSimulation` object.
"""
function read_Newtonian_KAWCA_simulation(filenamefull, dt; maxt=-1, every=1, original=false)
    println("Reading data file")
    f = open(filenamefull)
    r = Array{Array{Float64, 2}, 1}()
    v = Array{Array{Float64, 2}, 1}()
    types_array = Array{Array{Int, 1}, 1}()
    t = Vector{Float64}()
    iter = eachline(f)
    box_size = 0
    timestep = -1
    for line in iter
        if line == "ITEM: TIMESTEP"
            timestep += 1

            if length(r) >= maxt && maxt > 0
                break
            end
            push!(t, parse(Int64, iterate(iter)[1]))
            if timestep % every != 0
                continue
            end
            if length(r) % 100 == 0
                println(length(r), " snapshots found.")
            end
            _ = iterate(iter) # ITEM:NUMBER OF ATOMS
            N = parse(Int64, iterate(iter)[1])
            _ = iterate(iter) # ITEM: BOX BOUNDS
            box_size = parse(Float64, split(iterate(iter)[1])[2])
            _ = iterate(iter)
            _ = iterate(iter)
            _ = iterate(iter) # ITEM: ATOMS
            rnew = zeros(N, 3)
            vnew = zeros(N, 3)
	        types = zeros(Int, N)
            for i = 1:N
                index = Parsers.parse(Int64, readuntil(f, ' '))
                types[index] = Parsers.parse(Int64, readuntil(f, ' '))
                # if timestep == 0
                #     readuntil(f, ' ')
                # end

                for j = 1:3
                    rnew[index, j] = Parsers.parse(Float64, readuntil(f, ' '))
                end
                for j = 1:2
                    vnew[index, j] = Parsers.parse(Float64, readuntil(f, ' '))  
                end
                vnew[index, 3] = Parsers.parse(Float64, readline(f))  
                # readline(f)
            end
            push!(r, rnew)
            push!(v, vnew)
            if length(types_array) == 0
                 push!(types_array, types)
            end
        end
    end
    close(f)
    box_sizes = [box_size, box_size, box_size]
    t .-= 1
    println(length(r), " timesteps found with ", size(r[1])[1], " atoms.")
    r = reshape_data(r)
    remap_positions!(r, box_sizes)
    if original
        r = find_original_trajectories(r, box_sizes)
    end
    v = reshape_data(v)
    types = types_array[1]
    dt_arr = Int.([t[i]-t[1] for i in 1:length(t)])
    t1_t2_pair_array = [[1;; i] for i in 1:length(t)]
    F = zeros(size(v)...)
    rvec = separate_trajectories(r, types)
    vvec = separate_trajectories(v, types)
    Fvec = separate_trajectories(F, types)
    s = MultiComponentSimulation(sum(size.(rvec,2)), 
                            3, 
                            length(rvec), 
                            length(t), 
                            dt,
                            size.(rvec,2),
                            rvec,
                            vvec,
                            Fvec,
                            t*dt,
                            box_sizes,
                            dt_arr*dt,
                            t1_t2_pair_array,
                            filenamefull
                            )
    return s
end

"""
    read_Brownian_KALJ_simulation(filenamefull, dt; maxt=-1, every=1, original=false, forces=true)

Read a Brownian simulation of a Kob-Andersen Lennard-Jones mixture. The file was generated with LAMMPS. 
The columns in the data file correspond to the following:
1. Atom ID
2. Atom type
3. x position
4. y position
5. z position

# Arguments
- `filenamefull::String`: Path to the simulation file.
- `dt::Float64`: Time step.
- `maxt::Int=-1`: Maximum number of time steps to read.
- `every::Int=1`: Read every `every`-th time step.
- `original::Bool=false`: Whether to reconstruct the original trajectories.
- `forces::Bool=true`: Whether to compute forces.

# Returns
- `MultiComponentSimulation`: A `MultiComponentSimulation` object.
"""
function read_Brownian_KALJ_simulation(filenamefull, dt; maxt=-1, every=1, original=false, forces=true)
    println("Reading data file")
    f = open(filenamefull)
    r = Array{Array{Float64, 2}, 1}()
    # v = Array{Array{Float64, 2}, 1}()
    types_array = Array{Array{Int, 1}, 1}()
    t = Vector{Float64}()
    iter = eachline(f)
    box_size = 0
    timestep = -1
    for line in iter
        if line == "ITEM: TIMESTEP"
            timestep += 1

            if length(r) >= maxt && maxt > 0
                break
            end
            push!(t, parse(Int64, iterate(iter)[1]))
            if timestep % every != 0
                continue
            end
            if length(r) % 100 == 0
                println(length(r), " snapshots found.")
            end
            _ = iterate(iter) # ITEM:NUMBER OF ATOMS
            N = parse(Int64, iterate(iter)[1])
            _ = iterate(iter) # ITEM: BOX BOUNDS
            box_size = parse(Float64, split(iterate(iter)[1])[2])
            _ = iterate(iter)
            _ = iterate(iter)
            _ = iterate(iter) # ITEM: ATOMS
            rnew = zeros(N, 3)
            # vnew = zeros(N, 3)
	        types = zeros(Int, N)
            for i = 1:N
                index = Parsers.parse(Int64, readuntil(f, ' '))
                types[index] = Parsers.parse(Int64, readuntil(f, ' '))

                for j = 1:2
                    rnew[index, j] = Parsers.parse(Float64, readuntil(f, ' '))
                end
                rnew[index, 3] = Parsers.parse(Float64, readline(f))  


            end
            push!(r, rnew)
            if length(types_array) == 0
                 push!(types_array, types)
            end
        end
    end
    close(f)
    box_sizes = [box_size, box_size, box_size]
    println(length(r), " timesteps found with ", size(r[1])[1], " atoms.")
    r = reshape_data(r)
    remap_positions!(r, box_sizes)
    if original
        r = find_original_trajectories(r, box_sizes)
    end
    types = types_array[1]
    dt_arr, t1_t2_pair_array = find_allowed_t1_t2_pair_array_quasilog(t;doublefactor=10)
    #dt_arr = Int.([t[i]-t[1] for i in 1:length(t)])
    #t1_t2_pair_array = [[1;; i] for i in 1:length(t)]
    F = zeros(size(r)...)
    v = zeros(size(r)...)

    rvec = separate_trajectories(r, types)
    vvec = separate_trajectories(v, types)
    Fvec = separate_trajectories(F, types)
    U = KAWCA(1.0, 1.5,0.5,1.0,0.8,0.88)

    s = MultiComponentSimulation(sum(size.(rvec,2)), 
                            3, 
                            length(rvec), 
                            length(t),
                            dt,
                            size.(rvec,2),
                            rvec,
                            vvec,
                            Fvec,
                            t*dt,
                            box_sizes,
                            dt_arr*dt,
                            t1_t2_pair_array,
                            filenamefull
                            )
    if forces
        println("Computing Forces")
        calculate_forces!(s, U; cutoff=2.5)
    end
    return s
end


"""
    read_Brownian_KAWCA_simulation(filenamefull, dt; maxt=-1, every=1, original=false, forces=true)

Read a Brownian simulation of a Kob-Andersen WCA mixture. The file was generated with LAMMPS. 
The columns in the data file correspond to the following:
1. Atom ID
2. Atom type
3. x position
4. y position
5. z position

# Arguments
- `filenamefull::String`: Path to the simulation file.
- `dt::Float64`: Time step.
- `maxt::Int=-1`: Maximum number of time steps to read.
- `every::Int=1`: Read every `every`-th time step.
- `original::Bool=false`: Whether to reconstruct the original trajectories.
- `forces::Bool=true`: Whether to compute forces.

# Returns
- `MultiComponentSimulation`: A `MultiComponentSimulation` object.
"""
function read_Brownian_KAWCA_simulation(filenamefull, dt; maxt=-1, every=1, original=false, forces=true)
    println("Reading data file")
    f = open(filenamefull)
    r = Array{Array{Float64, 2}, 1}()
    # v = Array{Array{Float64, 2}, 1}()
    types_array = Array{Array{Int, 1}, 1}()
    t = Vector{Float64}()
    iter = eachline(f)
    box_size = 0
    timestep = -1
    for line in iter
        if line == "ITEM: TIMESTEP"
            timestep += 1

            if length(r) >= maxt && maxt > 0
                break
            end
            push!(t, parse(Int64, iterate(iter)[1]))
            if timestep % every != 0
                continue
            end
            if length(r) % 100 == 0
                println(length(r), " snapshots found.")
            end
            _ = iterate(iter) # ITEM:NUMBER OF ATOMS
            N = parse(Int64, iterate(iter)[1])
            _ = iterate(iter) # ITEM: BOX BOUNDS
            box_size = parse(Float64, split(iterate(iter)[1])[2])
            _ = iterate(iter)
            _ = iterate(iter)
            _ = iterate(iter) # ITEM: ATOMS
            rnew = zeros(N, 3)
            # vnew = zeros(N, 3)
	        types = zeros(Int, N)
            for i = 1:N
                index = Parsers.parse(Int64, readuntil(f, ' '))
                types[index] = Parsers.parse(Int64, readuntil(f, ' '))

                for j = 1:2
                    rnew[index, j] = Parsers.parse(Float64, readuntil(f, ' '))
                end
                rnew[index, 3] = Parsers.parse(Float64, readline(f))  


            end
            push!(r, rnew)
            if length(types_array) == 0
                 push!(types_array, types)
            end
        end
    end
    close(f)
    box_sizes = [box_size, box_size, box_size]
    println(length(r), " timesteps found with ", size(r[1])[1], " atoms.")
    r = reshape_data(r)
    remap_positions!(r, box_sizes)
    if original
        r = find_original_trajectories(r, box_sizes)
    end
    types = types_array[1]
    dt_arr, t1_t2_pair_array = find_allowed_t1_t2_pair_array_quasilog(t; doublefactor=10 )
    #dt_arr = Int.([t[i]-t[1] for i in 1:length(t)])
    #t1_t2_pair_array = [[1;; i] for i in 1:length(t)]
    F = zeros(size(r)...)
    v = zeros(size(r)...)

    rvec = separate_trajectories(r, types)
    vvec = separate_trajectories(v, types)
    Fvec = separate_trajectories(F, types)
    U = KAWCA(1.0, 1.5,0.5,1.0,0.8,0.88)

    s = MultiComponentSimulation(sum(size.(rvec,2)), 
                            3, 
                            length(rvec), 
                            length(t),
                            dt,
                            size.(rvec,2),
                            rvec,
                            vvec,
                            Fvec,
                            t*dt,
                            box_sizes,
                            dt_arr*dt,
                            t1_t2_pair_array,
                            filenamefull
                            )
    if forces
        println("Computing Forces")
        calculate_forces!(s, U; cutoff=2.0^(1.0/6.0))
    end
    return s
end

"""
    read_monodisperse_hard_sphere_simulation(filename; original=false, velocities=false, forcestype=false, dtarr=true)

Read a simulation of monodisperse particles. The file was generated with SimulationCode.jl.

# Arguments
- `filename::String`: Path to the simulation file.
- `original::Bool=false`: Whether to reconstruct the original trajectories.
- `velocities::Bool=false`: Whether to read velocities.
- `forcestype=false`: Type of forces to compute.
- `dtarr::Bool=true`: Whether to compute the time step array.

# Returns
- `SingleComponentSimulation`: A `SingleComponentSimulation` object.
"""
function read_monodisperse_hard_sphere_simulation(filename; original=false, velocities=false, forcestype=false, dtarr=true)
    # println("Reading data file")
    f = h5open(filename)
    saved_at_times = sort(parse.(Int64, keys(f["positions"])))
    Ndims, N = size(read(f["positions"][string(saved_at_times[1])]))

    
    box_size = read(HDF5.attributes(f)["box_size"])
    dt = float(read(HDF5.attributes(f)["Δt"]))
    dt = ifelse(dt == 0.0, 1.0, dt)
    m = 0.0
    try
        m += read(HDF5.attributes(f)["m"])
    catch
        m += 1.0
    end
    kBT = float(read(HDF5.attributes(f)["kBT"]))


    r = zeros(Ndims, N, length(saved_at_times))
    if velocities
         v = zeros(Ndims, N, length(saved_at_times))
    else
         v = zeros(1,1,1)
    end
    if forcestype != false
         F = zeros(Ndims, N, length(saved_at_times))
    else
         F = zeros(1,1,1)
    end
    D = zeros(N)

    if "diameters" in keys(f["diameters"])
        D .= read(f["diameters"]["diameters"])
    else
        D .= read(f["diameters"]["0"])
    end

    N_written = 0
    for t in saved_at_times
        N_written += 1
        r[:, :, N_written] .= read(f["positions"][string(t)])
        if velocities
            v[:, :, N_written] .= read(f["velocities"][string(t)])
        end
    end

    close(f)
    box_sizes = [box_size for i = 1:Ndims]
    if original
        r = find_original_trajectories(r, box_sizes)
    end
    if dtarr
        dt_arr, t1_t2_pair_array = find_allowed_t1_t2_pair_array_quasilog(saved_at_times; doublefactor=10)
    else
        dt_arr, t1_t2_pair_array = [1, 2], [zeros(Int, 2,2)]
    end
    s = SingleComponentSimulation(N, Ndims, length(saved_at_times), dt, r, v, F, D, saved_at_times*dt, box_sizes, dt_arr*dt, t1_t2_pair_array, filename)
    calculate_forces!(s, forcestype)
    return s
end

"""
    read_simulation_Berthier(filename; original=false, velocities=false, forcestype=false, time_origins="quasilog")

Read a simulation of the polydisperse system by Berthier. The file was generated with SimulationCode.jl.

# Arguments
- `filename::String`: Path to the simulation file.
- `original::Bool=false`: Whether to reconstruct the original trajectories.
- `velocities::Bool=false`: Whether to read velocities.
- `forcestype=false`: Type of forces to compute.
- `time_origins="quasilog"`: How to select the time origins.

# Returns
- `SingleComponentSimulation`: A `SingleComponentSimulation` object.
"""
function read_simulation_Berthier(filename; original=false, velocities=false, forcestype=false, time_origins="quasilog")
    println("Reading data file")
    f = h5open(filename)
    saved_at_times = sort(parse.(Int64, keys(f["positions"])))
    Ndims, N = size(read(f["positions"][string(saved_at_times[1])]))
    
    box_size = read(HDF5.attributes(f)["box_size"])
    dt = float(read(HDF5.attributes(f)["Δt"]))
    dt = ifelse(dt == 0, 1.0, dt) 
    m = 0.0
    try
        m += read(HDF5.attributes(f)["m"])
    catch
        m += 1.0
    end
    kBT = float(read(HDF5.attributes(f)["kBT"]))


    r = zeros(Ndims, N, length(saved_at_times))
    if velocities
         v = zeros(Ndims, N, length(saved_at_times))
    else
         v = zeros(1,1,1)
    end
    if forcestype != false
         F = zeros(Ndims, N, length(saved_at_times))
    else
         F = zeros(1,1,1)
    end
    D = zeros(N, length(saved_at_times))
    D = zeros(N)
    D[:] .= read(f["diameters"]["diameters"])
    N_written = 0
    for t in saved_at_times
        N_written += 1
        r[:, :, N_written] .= read(f["positions"][string(t)])
        if velocities
            v[:, :, N_written] .= read(f["velocities"][string(t)])
        end
    end

    close(f)
    box_sizes = [box_size for i = 1:Ndims]
    if original
        r = find_original_trajectories(r, box_sizes)
    end


    if time_origins == "quasilog"
        dt_arr, t1_t2_pair_array = find_allowed_t1_t2_pair_array_quasilog(saved_at_times; doublefactor=10)
    elseif typeof(time_origins) == Int
        dt_arr, t1_t2_pair_array =  find_allowed_t1_t2_pair_array_log_multstarts(saved_at_times, time_origins)
    else
        error("Specify time origins")
    end
    s = SingleComponentSimulation(N, Ndims, length(saved_at_times), dt, r, v, F, D, saved_at_times*dt, box_sizes, dt_arr*dt, t1_t2_pair_array, filename)
    calculate_forces!(s, forcestype)
    return s
end


"""
    reshape_data(r::Array{Array{Float64,2},1})

Reshape the data from a vector of 2D arrays to a 3D array.

# Arguments
- `r::Array{Array{Float64,2},1}`: A vector of 2D arrays.

# Returns
- `Array{Float64,3}`: A 3D array.
"""
function reshape_data(r::Array{Array{Float64,2},1})
    N_timesteps = length(r)
    N, Ndim = size(r[1])    
    rnew = zeros(Ndim, N, N_timesteps)
    for t = 1:N_timesteps
        for particle = 1:N
            for dim = 1:Ndim
                rnew[dim, particle, t] = r[t][particle, dim]
            end
        end
    end
    return rnew
end

"""
    separate_trajectories(r, type_list)

Separate the trajectories of different species.

# Arguments
- `r`: The trajectories.
- `type_list`: A list of particle types.

# Returns
- A vector of trajectories, one for each species.
"""
function separate_trajectories(r, type_list)
    N_species = length(unique(type_list))
    Ndim, N, N_timesteps = size(r)
    @assert minimum(type_list) == 1
    @assert maximum(type_list) == N_species

    N_particles_per_species = zeros(Int, N_species)
    for type in type_list
        N_particles_per_species[type] += 1
    end
    rvec = [zeros(Ndim, N_particles, N_timesteps) for N_particles in N_particles_per_species]
    particles_done = zeros(Int, N_species)
    for particle = 1:N
        species = type_list[particle]
        particles_done[species] += 1
        for timestep = 1:N_timesteps
            for dim = 1:Ndim
                rvec[species][dim, particles_done[species], timestep] = r[dim, particle, timestep]
            end
        end

    end
    return rvec
end

"""
    remap_positions!(r::Array{Float64,3}, box_sizes)

Remap the positions to be inside the box.

# Arguments
- `r::Array{Float64,3}`: The particle positions.
- `box_sizes`: The box sizes.
"""
function remap_positions!(r::Array{Float64,3}, box_sizes)
    N_timesteps = length(r)
    Ndim, N, N_timesteps = size(r)
    for t = 1:N_timesteps
        for particle = 1:N
            for dim = 1:Ndim
                box_size = box_sizes[dim]
                r[dim, particle, t] -= floor(r[dim, particle, t]/box_size)*box_size
            end
        end
    end
end


"""
    find_original_trajectories(r, box_sizes)

Reconstruct the original trajectories from the remapped ones.

# Arguments
- `r`: The remapped trajectories.
- `box_sizes`: The box sizes.

# Returns
- The original trajectories.
"""
function find_original_trajectories(r, box_sizes)
    Ndim, N, N_timesteps = size(r)  
    rr = deepcopy(r)
    for particle = 1:N
        for dim = 1:Ndim
            box_size = box_sizes[dim]
            for t = 1:N_timesteps-1
                if rr[dim, particle, t+1] - rr[dim, particle, t] > box_size/2
                    for t2 = t+1:N_timesteps
                        rr[dim, particle, t2] = rr[dim, particle, t2] - box_size
                    end
                elseif rr[dim, particle, t+1] - rr[dim, particle, t] < -box_size/2
                    for t2 = t+1:N_timesteps
                        rr[dim, particle, t2] =rr[dim, particle, t2] + box_size
                    end                
                end
            end
        end
    end
    return rr
end

"""
    find_quasilog_time_array(maxsteps; doublefactor=10)

Create a quasi-logarithmic time array.
eg if double_factor is 10, the dt array will be 
    [1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,200,...]

# Arguments
- `maxsteps`: The maximum number of steps.
- `doublefactor=10`: The factor by which to double the time step.

# Returns
- A quasi-logarithmic time array.
"""
function find_quasilog_time_array(maxsteps; doublefactor=10)
    save_array = Int64[]
    t = 0
    dt = 1
    while t <= maxsteps
        if !(t % (10*dt) == 0) || t == 0 
            push!(save_array, t)
        end
        t += dt
        if t == dt*doublefactor
            dt *= 10
        end
    end
    return save_array
end

"""
    find_allowed_t1_t2_pair_array_quasilog(t_array; doublefactor=150, dt=1.0)

Find allowed time pairs for a quasi-logarithmic time array.


# Arguments
- `t_array`: The time array.
- `doublefactor=150`: The factor by which to double the time step.

# Returns
- A tuple containing the `dt_array` and `t1_t2_pair_array`.
"""
function find_allowed_t1_t2_pair_array_quasilog(t_array::Vector{Int}; doublefactor=150)
    maxt = t_array[end]-t_array[1]
    dt_array = find_quasilog_time_array(maxt; doublefactor=doublefactor)
    t1_t2_pair_array = Vector{Array{Int64, 2}}()
    for dt in dt_array
        tstart_arr = zeros(Int64, 0, 2)
        for (t1idx, tstart) in enumerate(t_array)
            t2indx = findfirst(isequal(tstart+dt), t_array)
            if !isnothing(t2indx)
                tstart_arr = cat(tstart_arr, [t1idx t2indx], dims=1)
            end
        end
        push!(t1_t2_pair_array, tstart_arr)
    end
    return dt_array, t1_t2_pair_array
end

"""
    find_log_time_array_multiple_starts(log_factor, N_starts, N_max)

Create a logarithmic time array with multiple starts.

# Arguments
- `log_factor`: The logarithmic factor.
- `N_starts`: The number of starting points.
- `N_max`: The maximum number of steps.

# Returns
- A logarithmic time array with multiple starts.
"""
function find_log_time_array_multiple_starts(log_factor, N_starts, N_max)
    start_times = 0:(N_max÷N_starts):N_max
    when_to_save = Int[collect(start_times)...]
    for i_start in start_times
        t = 1
        while t <= N_max
            push!(when_to_save, t+i_start)
            t *= log_factor
            t = ceil(Int, t)
        end
    end
    push!(when_to_save, N_max)
    return sort(unique(when_to_save[when_to_save .<= N_max]))
end


"""
    find_allowed_t1_t2_pair_array_log_multstarts(t_array, N_starts)

Find allowed time pairs for a logarithmic time array with multiple starts.

# Arguments
- `t_array`: The time array in units of dt.
- `N_starts`: The number of starting points.

# Returns
- A tuple containing the `dt_array` and `t1_t2_pair_array`.
"""
function find_allowed_t1_t2_pair_array_log_multstarts(t_array::Vector{Int}, N_starts)
    maxt = t_array[end]

    dt_array = find_log_time_array_multiple_starts(1.3, 1, maxt)
    @assert all(dt in t_array for dt in dt_array)

    t1_t2_pair_array = Vector{Array{Int64, 2}}()
    for dt in dt_array
        tstart_arr = Vector{Vector{Int64}}()
        for t1 in 0:(maxt÷N_starts):maxt
            t2 = t1 + dt
            if t2 > maxt
                break
            end
            @assert t2 in t_array
            @assert t1 in t_array
            it1 = findfirst(isequal(t1), t_array)
            it2 = findfirst(isequal(t2), t_array)
            push!(tstart_arr, [it1, it2])
        end
        push!(t1_t2_pair_array, stack(tstart_arr)')
    end
    return dt_array, t1_t2_pair_array
end



"""
    Find allowed time pairs for any time array.

Loops through all pairs of times in the array and finds the time differences that are present in the array.

# Arguments
- `t_array`: The time array.

# Returns
- A vector of allowed time differences.
- A vector of allowed time pairs for every dt.


"""
function find_allowed_dt_array(t_array::Vector{Int})
    allowed_dt_array = Vector{Int}()

    for t1 in t_array
        for t2 in t_array
            if t2 > t1
                dt = t2 - t1
                if !(dt in allowed_dt_array)
                    push!(allowed_dt_array, dt)
                end
            end
        end
    end

    dt_array = sort(collect(unique(allowed_dt_array)))
    t1_t2_pair_array = Vector{Array{Int64, 2}}()
    for dt in dt_array
        t1_t2s = Vector{Vector{Int64}}()
        for t1 in t_array
            t1_index = findfirst(isequal(t1), t_array)
            t2 = t1 + dt
            t2_index = findfirst(isequal(t2), t_array)
            if t2 in t_array
                push!(t1_t2s, [t1_index, t2_index])
            end
        end
        push!(t1_t2_pair_array, stack(t1_t2s)')
    end
    return dt_array, t1_t2_pair_array
end



"""
    calculate_forces!(s, forcestype::Bool)

Calculate the forces for a simulation.

This is a placeholder function that throws an error if `forcestype` is true.

# Arguments
- `s`: The simulation object.
- `forcestype::Bool`: Whether to calculate forces.
"""
calculate_forces!(s, forcestype::Bool) = forcestype ? error("Specify force type") : nothing