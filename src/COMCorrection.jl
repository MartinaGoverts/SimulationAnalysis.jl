"""
    COM_correction_function(r,box_sizes, N, original::Bool=false)

Updates the positions by applying a COM correction.

# Arguments
    - `r`: A three dimensional array containing all saved particle positions. Size is 
    given by `d x N x Nt`, where `N` is the number of particles, and `Nt` the number of frames
    - `box_sizes`: A tuple containing the sizes of the simulation box
    - `original::Bool=false`: Whether to reconstruct the original trajectories.

# Returns
    - `updated_r_array`: An array containing all updated positions.
    - if original = true: 
    `updated_r_array_original`: An array containing all updated positions with reconstructed original trajectories.
"""


function COM_correction_function(r, box_sizes, original)
    # first apply COM correction on original, then convert back to normal by applying pbc?
    # required to avoid issues where particles end up outside the box!
    N = size(r, 2)
    r_original = find_original_trajectories(r, box_sizes)

    # find COM
    first_COM = sum(r_original[:, :, 1], dims=2) ./ N

    updated_r_array_original = zeros(size(r_original))

    # Update positions
    for i in 1:size(r_original,3)
        shift = sum(r_original[:,:,i], dims=2) ./ N .- first_COM
        updated_r_array_original[:,:,i] .= r_original[:,:,i] .- shift
    end

    if original
        return updated_r_array_original
    else
    return apply_periodic_boundary_conditions(updated_r_array_original, box_sizes)
    end
end

"""
apply_periodic_boundary_conditions(position, box_sizes)

Applies periodic boundary conditions to a given position, ensuring it wraps around the simulation box.
For a position coordinate `x` and a box dimension `L`, the new coordinate `x'` is `x - floor(x/L) * L`.
This maps `x` to the interval `[0, L)`. The same logic applies to all dimensions.

# Arguments
- `position`: The original position, typically an `SVector` or any `AbstractVector` representing coordinates (e.g., `[x, y]`). It's first dimension must be equal to the length of `box_sizes`.
- `box_sizes`: The dimensions of the simulation box, typically an `SVector` or `AbstractVector` (e.g., `[Lx, Ly]`).

# Returns
- `new_position`: The position after applying periodic boundary conditions, of the same type as `position`.
"""

function apply_periodic_boundary_conditions(position, box_sizes)
    new_position = position .- floor.(position ./ box_sizes) .* box_sizes
    return new_position
end