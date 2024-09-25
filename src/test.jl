function find_k_array(r_array, box_size, k_max)
    N_dims = size(r_array, 1)
    dk = 2π/box_size
    N_k = round(Int, k_max/dk)
    k_array = zeros(N_dims, (N_k+1)^N_dims)
    ik = 1
    # these loops only work for 3d:
    for ikx in 0:N_k
        for iky in 0:N_k
            for ikz in 0:N_k
                k_array[1, ik] = ikx*dk
                k_array[2, ik] = iky*dk
                k_array[3, ik] = ikz*dk
                ik += 1
            end
        end
    end
    k_lengths = sqrt.(sum(k_array.^2, dims=1))[:]
    return k_array, k_lengths
end

function find_structure_factor(r_array, k_array)
    N_k = size(k_array, 2)
    N_dims, N_particles, N_frames = size(r_array)
    ρₖ = zeros(Complex{Float64}, N_k, N_frames) # matrix of complex numbers
    println("Calculating the structure factor for $N_particles particles, at $N_frames frames, for $N_k wave vectors.")
    for frame in 1:N_frames
        @show frame
        for ik in 1:N_k
            for particle in 1:N_particles
                kr = 0.0
                for dim in 1:N_dims
                    kr += r_array[dim, particle, frame]*k_array[dim, ik]
                end
                ρₖ[ik, frame] += exp(1im*kr)
            end
        end
    end
    S = sum(conj.(ρₖ) .* ρₖ, dims=2)/N_frames/N_particles
    S[1] = 0.0
    return real.(S)
end

findnearest(arr,val) = findmin(abs.(arr .- val))[2]

function bin_structure_factor(k_lengths, S, k_bin_array)
    N_bins = length(k_bin_array)
    S_binned = zeros(N_bins)
    count = zeros(N_bins)
    for (ik,k) in enumerate(k_lengths)
        ik_sample = findnearest(k_bin_array, k)
        S_binned[ik_sample] += S[ik]
        count[ik_sample] += 1
    end
    return S_binned ./ count
end

# generate random particle positions
N_dims = 3
N_particles = 1000
N_frames = 10
box_size = 10.0
r_array = rand(N_dims, N_particles, N_frames)*box_size

# find the structure factor
k_max = 15.0
k_array, k_lengths = find_k_array(r_array, box_size, k_max)
S = find_structure_factor(r_array, k_array)

# bin the structure factor
k_bin_array = range(0, k_max, length=30)
S = bin_structure_factor(k_lengths, S, k_bin_array)

## to install the plotting package run:
# import Pkg
# Pkg.add("Plots")

# plot the structure factor
using Plots
scatter(k_sample_array, S, xlabel="k", ylabel="S(k)", legend=false)