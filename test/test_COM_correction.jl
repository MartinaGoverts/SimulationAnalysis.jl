function test_shifts(r_original, N)
    # COM in X direction
    first_x_positions = r_original[1,:,1]
    first_COM_x = sum(first_x_positions)/N

    first_y_positions = r_original[2,:,1]
    first_COM_y = sum(first_y_positions)/N

    x_shift = []
    y_shift = []

    for i in range(1,size(r_original)[3])
        # x correction
        x_position = r_original[1,:,i]
        x_shift = push!(x_shift,sum(x_position)/N - first_COM_x)

        # y correction
        y_position = r_original[2,:,i]
        y_shift = push!(y_shift,sum(y_position)/N - first_COM_y)
    end
    return x_shift, y_shift
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
r_test = SimulationAnalysis.COM_correction_function(r_array, box_sizes, N, true)
x,y = test_shifts(r_test, N)
max_val = max(maximum(x), maximum(y))
@test max_val < 1e-10

# Check if all particles are still in the box by applying the COM
# correction with original = false
r_test_2 = SimulationAnalysis.COM_correction_function(r_array, box_sizes, N, false) 
@test !any(>(box_sizes[1]), r_test_2[1,:,:])
@test !any(>(box_sizes[2]), r_test_2[2,:,:])