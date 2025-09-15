
# load the data (this requires the SelfPropelledVoronoi package)
traj, params = SelfPropelledVoronoi.load_trajectory("data/spv_data_2.h5");
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
