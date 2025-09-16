
# load the data (this requires the SelfPropelledVoronoi package) 
# traj, params = SelfPropelledVoronoi.load_trajectory("data/spv_data_2.h5");

# we cannot depend on the unregistered SelfPropelledVoronoi package in this package
# so we re-implement the load_trajectory function here
using HDF5
function load_trajectory(filename::String)
    h5open(filename, "r") do file
        # Read simulation parameters (similar to load_simulation_state)
        params_group = file["parameters"]
        N = read(params_group["N"])
        dt = read(params_group["dt"])
        kBT = read(params_group["kBT"])
        frictionconstant = read(params_group["frictionconstant"])
        
        box_sizes_vec = read(params_group["box_sizes"])
        box = (; box_sizes=(box_sizes_vec[1], box_sizes_vec[2]))

        particles_group = params_group["particles"]
        target_perimeters = read(particles_group["target_perimeters"])
        target_areas = read(particles_group["target_areas"])
        K_P = read(particles_group["K_P"])
        K_A = read(particles_group["K_A"])
        active_force_strengths = read(particles_group["active_force_strengths"])
        rotational_diffusion_constants = read(particles_group["rotational_diffusion_constants"])

        particles = (; target_perimeters=target_perimeters, target_areas=target_areas, K_P=K_P, K_A=K_A, active_force_strengths=active_force_strengths, rotational_diffusion_constants=rotational_diffusion_constants)

        dump_info = (; filename=filename) # Using defaults
        rng = Random.MersenneTwister()
        if "seed" in keys(params_group)
            seed = read(params_group["seed"])
            Random.seed!(rng, seed)
        end
        callback = nothing

        parameter_struct = (;
            N=N, dt=dt, kBT=kBT, frictionconstant=frictionconstant, 
            periodic_boundary_layer_depth=3.0, verbose=false, box=box,
            particles=particles, dump_info=dump_info, callback=callback, RNG=rng
        )

        # Initialize TrajectoryData
        trajectory_data = (; positions_trajectory = Vector{Vector{SVector{2, Float64}}}(), 
                             orientations_trajectory = Vector{Vector{Float64}}(), 
                             forces_trajectory = Vector{Vector{SVector{2, Float64}}}(), 
                             potential_energy_trajectory = Float64[], 
                             areas_trajectory = Vector{Vector{Float64}}(), 
                             perimeters_trajectory = Vector{Vector{Float64}}(), 
                             steps_saved = Int[])

        # Identify and sort step groups
        step_group_names = String[]
        for name in keys(file)
            if all(isdigit, name)
                push!(step_group_names, name)
            end
        end
        
        # Sort step groups numerically
        sort!(step_group_names, by=x->parse(Int, x))

        if isempty(step_group_names)
            @warn "No simulation step groups found in HDF5 file: $filename. Returning empty TrajectoryData."
            return trajectory_data, parameter_struct # Return empty trajectory and params
        end

        for step_str in step_group_names
            step_group = file[step_str]
            current_step = parse(Int, step_str)

            # Load and append positions
            positions_raw = read(step_group["positions"])
            loaded_positions = [SVector{2, Float64}(positions_raw[:, i]) for i in 1:N]
            push!(trajectory_data.positions_trajectory, loaded_positions)

            # Load and append orientations
            if "orientations" in keys(step_group)
                loaded_orientations = read(step_group["orientations"])
                push!(trajectory_data.orientations_trajectory, loaded_orientations)
            else
                # If not saved, append a vector of zeros (or handle as per desired behavior)
                push!(trajectory_data.orientations_trajectory, zeros(Float64, N))
            end

            # Load and append forces
            if "forces" in keys(step_group)
                forces_raw = read(step_group["forces"])
                loaded_forces = [SVector{2, Float64}(forces_raw[:,i]) for i in 1:N]
                push!(trajectory_data.forces_trajectory, loaded_forces)
            else
                # If not saved, append a vector of zero vectors
                push!(trajectory_data.forces_trajectory, [zeros(SVector{2, Float64}) for _ in 1:N])
            end

            # Load and append potential energy
            if "potential_energy" in keys(step_group)
                loaded_potential_energy = read(step_group["potential_energy"])
                push!(trajectory_data.potential_energy_trajectory, loaded_potential_energy)
            else
                push!(trajectory_data.potential_energy_trajectory, 0.0) # Default if not saved
            end
            
            # Load and append areas
            loaded_areas = read(step_group["areas"])
            push!(trajectory_data.areas_trajectory, loaded_areas)

            # Load and append perimeters
            loaded_perimeters = read(step_group["perimeters"])
            push!(trajectory_data.perimeters_trajectory, loaded_perimeters)
            
            # Append the current step number
            push!(trajectory_data.steps_saved, current_step)
        end

        return trajectory_data, parameter_struct
    end
end

traj, params = load_trajectory("data/spv_data_2.h5");
sim = SimulationAnalysis.read_SPV_simulation(traj, params, original=true);

kbounds = (0, 20.0);
kspace = SimulationAnalysis.construct_k_space(sim, kbounds);
current_modes = SimulationAnalysis.find_current_modes(sim, kspace; verbose=false);

@test size(current_modes.Re) == size(current_modes.Im) == (sim.Nt, kspace.Nk)
@test current_modes.Re[1,1] == 5.863837291322234
@test current_modes.Im[100,1] == 4.349061145734272

# calculation of w(k)
@test SimulationAnalysis.find_static_velocity_correlations(sim, kspace, current_modes, kmin=0.0, kmax=20.0) ==
    SimulationAnalysis.real_static_correlation_function(current_modes.Re, current_modes.Im, current_modes.Re, current_modes.Im, kspace, 0.0, 20.0) / sim.N

k_array = 0.4:0.5:10;
wk_array = SimulationAnalysis.find_static_velocity_correlations(sim, kspace, current_modes, k_array)
@test any(isnan.(wk_array)) == false
@test maximum(abs.(wk_array)) <= 1.0

# set active forces to zero: total force == interaction force
sim2 = SelfPropelledVoronoiSimulation(
    sim.N,
    sim.Ndims,
    sim.Nt,
    sim.dt,
    0.0,
    sim.mobility,
    sim.r_array,
    sim.u_array,
    sim.F_array,
    sim.perimeter_array,
    sim.area_array,
    sim.Epot_array,
    sim.t_array,
    sim.box_sizes,
    sim.dt_array,
    sim.t1_t2_pair_array,
    ""
)

@test sim2.v0 == 0.0
@test sim2.F_array == SimulationAnalysis.calculate_total_force(sim2)

# set interaction force to zero: total force == active force
sim3 = SelfPropelledVoronoiSimulation(
    sim.N,
    sim.Ndims,
    sim.Nt,
    sim.dt,
    sim.v0,
    sim.mobility,
    sim.r_array,
    sim.u_array,
    zero(sim.F_array),
    sim.perimeter_array,
    sim.area_array,
    sim.Epot_array,
    sim.t_array,
    sim.box_sizes,
    sim.dt_array,
    sim.t1_t2_pair_array,
    ""
)

Ftot3 = SimulationAnalysis.calculate_total_force(sim3);
@test sim3.F_array != Ftot3
@test Ftot3[1,:,:] == sim3.v0 / sim3.mobility * cos.(sim.u_array)
