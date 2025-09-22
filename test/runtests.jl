using Test, SimulationAnalysis, Random, StaticArrays

file = joinpath(@__DIR__, "data", "test_trajectory.h5")

# Load the trajectory



# tests

@testset "read trajectory" begin
    traj = SimulationAnalysis.read_simulation_Berthier(file; original=false, velocities=false, forcestype=false, time_origins=10)

    @test traj.N == 1000
    @test size(traj.r_array) == (3, 1000, 430) 
    @test size(traj.v_array) == (1, 1, 1)
    @test size(traj.F_array) == (1, 1, 1)

    @test size(traj.D_array) == (1000, )
    @test size(traj.t_array) == (430, )
end

@testset "KSpace" begin
    traj = SimulationAnalysis.read_simulation_Berthier(file; original=false, velocities=false, forcestype=false, time_origins=10)

    kspace = SimulationAnalysis.construct_k_space(traj, (0.0, 3.0); kfactor=1, negative=true, rectangular=false)
    kspace = SimulationAnalysis.construct_k_space(traj, (0.0, 3.0); kfactor=1, negative=true, rectangular=true)
    kspace = SimulationAnalysis.construct_k_space(traj, (0.0, 3.0); kfactor=1, negative=false, rectangular=false)

    # compute structure factor

    S = SimulationAnalysis.find_structure_factor(traj; kmin=2.0, kmax=2.4, kfactor=1)

    # compute ISF

    ISF = SimulationAnalysis.find_intermediate_scattering_function(traj; kmin=2.0, kmax=2.4, kfactor=1)

    @test abs(ISF[1] - S) < 1e-2

    # compute self ISF

    sISF = SimulationAnalysis.find_self_intermediate_scattering_function(traj, kspace; kmin=2.0, kmax=2.4)

end


@testset "g(r)" begin
    traj = SimulationAnalysis.read_simulation_Berthier(file; original=false, velocities=false, forcestype=false, time_origins=10)

    gr =  SimulationAnalysis.find_radial_distribution_function(traj, 10,  10.0)
    @test length(gr[1]) == length(gr[2])
    @test all(gr[2] .>= 0)
end

@testset "Neighborlists" begin
    traj = SimulationAnalysis.read_simulation_Berthier(file; original=false, velocities=false, forcestype=false, time_origins=10)

    neighborlists = @time SimulationAnalysis.find_absolute_distance_neighborlists(traj, 1.2)

    @test length(neighborlists) == 430

    @test length(neighborlists[1]) == 1000


    neighborlists = @time SimulationAnalysis.find_voronoi_neighborlists(traj; indices=1:10)

    @test length(neighborlists) == 430

    @test length(neighborlists[1]) == 1000


    # CB

    Cb = SimulationAnalysis.find_CB(traj, neighborlists, neighborlists)

    @test size(Cb) == (47, 1000)

end


@testset "MSD" begin
    traj = SimulationAnalysis.read_simulation_Berthier(file; original=false, velocities=false, forcestype=false, time_origins=10)

    MSD = SimulationAnalysis.find_mean_squared_displacement(traj)

    @test length(MSD) == length(traj.dt_array)
end

@testset "velocity correlations" begin
    include("test_velocity_correlations.jl")
end

@testset "COM correction" begin
    include("test_COM_correction.jl")
end