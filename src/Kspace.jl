
"""
    KSpace{dims, S, OA}

A struct to hold the k-space information of a simulation.

# Fields
- `s::S`: The simulation object.
- `Nk::Int64`: The number of k-vectors.
- `k_lengths::Array{Float64, 1}`: The lengths of the k-vectors.
- `k_array::Array{Float64, 2}`: The k-vectors.
- `kfactor::Int64`: The k-factor.
- `cartesian_to_linear::OA`: A mapping from Cartesian to linear indices.
"""
struct KSpace{dims, S, OA}
    s::S
    Nk::Int64
    k_lengths::Array{Float64, 1}
    k_array::Array{Float64, 2}
    kfactor::Int64
    cartesian_to_linear::OA
end

"""
    construct_k_space(s::Simulation, bounds; kfactor=1, negative=false, rectangular=false)

Construct the k-space of a simulation.

# Arguments
- `s::Simulation`: The simulation object.
- `bounds`: The bounds of the k-space.
- `kfactor=1`: The k-factor.
- `negative=false`: Whether to include negative k-vectors.
- `rectangular=false`: Whether to use a rectangular k-space.

# Returns
- A `KSpace` object.
"""
function construct_k_space(s::Simulation, bounds; kfactor=1, negative=false, rectangular=false)
    k_lengths, k_array, cartesian_to_linear = find_k_array(bounds, s.box_sizes, s.Ndims; kfactor=kfactor, negative=negative, rectangular=rectangular)
    return KSpace{s.Ndims, typeof(s), typeof(cartesian_to_linear)}(s, length(k_lengths), k_lengths, k_array, kfactor, cartesian_to_linear)
end

"""
    find_k_array(bounds, box_sizes, dims; kfactor=1, rectangular=false, negative=false)

Find the k-vectors in a given k-space.

# Arguments
- `bounds`: The bounds of the k-space.
- `box_sizes`: The box sizes.
- `dims`: The number of dimensions.
- `kfactor=1`: The k-factor.
- `rectangular=false`: Whether to use a rectangular k-space.
- `negative=false`: Whether to include negative k-vectors.

# Returns
- A tuple containing the `k_lengths`, `k_array`, and `cartesian_to_linear`.
"""
function find_k_array(bounds, box_sizes, dims; kfactor=1, rectangular=false, negative=false)
    if dims==3
        return find_k_array_3D(bounds, box_sizes; kfactor=kfactor, rectangular=rectangular, negative=negative)
    elseif dims == 2
        return find_k_array_2D(bounds, box_sizes; kfactor=kfactor, rectangular=rectangular, negative=negative)
    else
        error("This is not implemented")
    end
end

"""
    find_k_array_3D(bounds, box_sizes; kfactor=1, rectangular=false, negative=false)

Find the k-vectors in a 3D k-space.

# Arguments
- `bounds`: The bounds of the k-space.
- `box_sizes`: The box sizes.
- `kfactor=1`: The k-factor.
- `rectangular=false`: Whether to use a rectangular k-space.
- `negative=false`: Whether to include negative k-vectors.

# Returns
- A tuple containing the `k_lengths`, `k_array`, and `cartesian_to_linear`.
"""
function find_k_array_3D(bounds, box_sizes; kfactor=1, rectangular=false, negative=false)
    #bounds
    @assert bounds[1] < bounds[2]
    if rectangular
        @assert bounds[1] == 0.0
    end

    kmin = bounds[1]
    kmax = bounds[2]
    kmin2 = kmin^2
    kmax2 = kmax^2
    # resolution
    dkx = 2π/box_sizes[1]*kfactor
    dky = 2π/box_sizes[2]*kfactor
    dkz = 2π/box_sizes[3]*kfactor

    # max number of kpoints in 1D
    Nmaxx = ceil(Int64, kmax/dkx)
    Nmaxy = ceil(Int64, kmax/dky)
    Nmaxz = ceil(Int64, kmax/dkz)

    Nminx = ifelse(negative, -Nmaxx, 0)
    Nminy = ifelse(negative, -Nmaxy, 0)
    Nminz = ifelse(rectangular, -Nmaxz, 0) # we can always take kz to be >= 0, to save a factor 2 in computations. if rectangular, we do take all

    #get all k values
    k_array = Array{Tuple{Float64, Float64, Float64}, 1}()
    cartesian_to_linear = OffsetArray(zeros(Int, Nmaxx-Nminx+1, Nmaxy-Nminy+1, Nmaxz-Nminz+1), Nminx:Nmaxx, Nminy:Nmaxy, Nminz:Nmaxz)
    ik = 0
    for iz = Nminz:Nmaxz
        for iy = Nminy:Nmaxy
            for ix = Nminx:Nmaxx
                if rectangular || kmin2 <  (dkx^2 *ix^2 + dky^2 *iy^2 + dkz^2 *iz^2) < kmax2
                    ik += 1
                    push!(k_array, (dkx*ix, dky*iy, dkz*iz))
                    cartesian_to_linear[ix, iy, iz] = ik
                end
            end
        end
    end
    k_array = Array(reinterpret(reshape, Float64, k_array))
    # calculate lengths
    k_lengths = sqrt.(reshape(sum(k_array.^2; dims=1), size(k_array)[2]))
    return k_lengths, k_array, cartesian_to_linear
end

"""
    find_k_array_2D(bounds, box_sizes; kfactor=1, rectangular=false, negative=false)

Find the k-vectors in a 2D k-space.

# Arguments
- `bounds`: The bounds of the k-space.
- `box_sizes`: The box sizes.
- `kfactor=1`: The k-factor.
- `rectangular=false`: Whether to use a rectangular k-space.
- `negative=false`: Whether to include negative k-vectors.

# Returns
- A tuple containing the `k_lengths`, `k_array`, and `cartesian_to_linear`.
"""
function find_k_array_2D(bounds, box_sizes; kfactor=1, rectangular=false, negative=false)
    #bounds
    @assert bounds[1] < bounds[2]
    if rectangular
        @assert bounds[1] == 0.0
    end

    kmin = bounds[1]
    kmax = bounds[2]
    kmin2 = kmin^2
    kmax2 = kmax^2
    # resolution
    dkx = 2π/box_sizes[1]*kfactor
    dky = 2π/box_sizes[2]*kfactor

    # max number of kpoints in 1D
    Nmaxx = ceil(Int64, kmax/dkx)
    Nmaxy = ceil(Int64, kmax/dky)

    Nminx = ifelse(negative, -Nmaxx, 0)
    Nminy = ifelse(rectangular, -Nmaxy, 0)

    #get all k values
    k_array = Array{Tuple{Float64, Float64}, 1}()
    cartesian_to_linear = OffsetArray(zeros(Int, Nmaxx-Nminx+1, Nmaxy-Nminy+1), Nminx:Nmaxx, Nminy:Nmaxy)
    ik = 0
    for iy = Nminy:Nmaxy
        for ix = Nminx:Nmaxx
            if rectangular || kmin2 <  (dkx^2 *ix^2 + dky^2 *iy^2) < kmax2
                ik += 1
                push!(k_array, (dkx*ix, dky*iy))
                cartesian_to_linear[ix, iy] = ik
            end
        end
    end
    k_array = Array(reinterpret(reshape, Float64, k_array))
    # calculate lengths
    k_lengths = sqrt.(reshape(sum(k_array.^2; dims=1), size(k_array)[2]))
    return k_lengths, k_array, cartesian_to_linear
end

# function construct_k_space_sampled(s::Simulation, k_sample_array::Vector{Float64}; k_binwidth=0.1, max_samples=typemax(Int))
#     k_lengths, k_array = find_k_array_sampled(s.box_sizes, k_sample_array, k_binwidth, max_samples)
#     return Kspace{s.Ndims, typeof(s)}(length(k_lengths), k_lengths, k_array, 1)
# end

# function find_k_array_sampled(box_sizes, k_sample_array, k_binwidth, max_samples)
#     k_array = zeros(length(box_sizes), 0)
#     k_lengths = zeros(1, 0)
#     for sample in k_sample_array
#         kbounds = (sample-k_binwidth/2, sample+k_binwidth/2)
#         k_lengths_new, k_array_new = find_k_array(kbounds, box_sizes; factor=1, rectangular=false, negative=true)
#         Nk_new = length(k_lengths_new)
#         if max_samples < Nk_new
#             randperm = randperm(Nk_new)[1:max_samples]
#             k_lengths_new = k_lengths_new[randperm]
#             k_array_new = k_array_new[3, randperm]
#         end
#         k_lengths = [klengths; k_lengths_new]
#         k_array_new = [klengths k_lengths_new]
#     end
#     return k_lengths, k_array
# end


