function test_shifts(r_original)
    # COM in X direction
    N = size(r_original, 2)
    first_COM = sum(r_original[:, :, 1], dims=2) ./ N

    shifts = []

    for i in range(1,size(r_original)[3])

        shift = sum(r_original[:,:,i], dims=2) ./ N .- first_COM

        push!(shifts,shift)

    end
    return shifts
end

N = 400
box_sizes = (20,20)
Nsteps = 10

# Generate random positions
r_array = Array{Float64}(undef, 2, N, Nsteps)
for i in 1:Nsteps
    # Each position is an SVector{2}, multiply by box size
    positions = [rand(SVector{2, Float64}) .* box_sizes for _ in 1:N]
    
    # Store positions into r_array[:, :, i]
    # Convert to 2 x N matrix for assignment
    r_array[:, :, i] = reduce(hcat, positions)
end

# Check if COM is correctly applied by calculating the COM correction first
# with original = true, and then calculating if the shifts are indeed 0 if we 
# apply a COM correction on the ouput again
r_test = SimulationAnalysis.COM_correction_function(r_array, box_sizes, true)
shifts = test_shifts(r_test)
max_val = maximum(map(maximum, shifts))
@test max_val < 1e-10

# Check if all particles are still in the box by applying the COM
# correction with original = false
r_test_2 = SimulationAnalysis.COM_correction_function(r_array, box_sizes, false) 
@test !any(>(box_sizes[1]), r_test_2[1,:,:])
@test !any(>(box_sizes[2]), r_test_2[2,:,:])