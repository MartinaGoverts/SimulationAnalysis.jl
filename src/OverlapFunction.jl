"""
    find_overlap_function(s::Union{SingleComponentSimulation, SelfPropelledVoronoiSimulation}; a=0.5)

Calculates the overlap function `Q(t)`.

The overlap function is a measure of the similarity between the particle configurations at time `0` and time `t`. It is defined as:
`Q(t) = (1/N) * Σ_i θ(a - |r_i(t) - r_i(0)|)`,
where `θ` is the Heaviside step function, `a` is a cutoff distance, and the average is taken over all particles `i` and time origins.
A value of 1 means the configurations are identical (within the cutoff), and a value of 0 means they are completely different.

# Arguments
- `s::Union{SingleComponentSimulation, SelfPropelledVoronoiSimulation}`: The simulation data.
- `a::Float64=0.5`: The cutoff distance for calculating the overlap. Typically this is a fraction of the particle diameter.

# Returns
- `Fs::Vector{Float64}`: A vector containing `Q(t)` for each time delay `Δt` in `s.dt_array`.
- `Fs_pp::Matrix{Float64}`: A `(Ndt, N)` matrix containing the overlap function for each particle.
"""
function find_overlap_function(s::Union{SingleComponentSimulation, SelfPropelledVoronoiSimulation}; a=0.5)
    N = s.N
    Ndt = length(s.dt_array)
    box_sizes = s.box_sizes
    dims = s.Ndims
    Fs = zeros(Ndt)
    Fs_pp = zeros(Ndt, N)
    for iδt in eachindex(s.dt_array)
        pairs_idt = s.t1_t2_pair_array[iδt]
        Npairs = size(pairs_idt, 1)
        for ipair = 1:Npairs
            t1 = pairs_idt[ipair, 1]
            t2 = pairs_idt[ipair, 2]
            for particle = 1:N
                dr2 = 0.0
                for dim = 1:dims
                    drdim = (s.r_array[dim, particle, t1] - s.r_array[dim, particle, t2])
                    drdim -= box_sizes[dim]*round(drdim/box_sizes[dim])
                    dr2 += drdim^2
                end
                r12 = sqrt(dr2)
                if r12 < a
                    Fs_pp[iδt, particle] += 1
                end
            end
        end
        Fs_pp[iδt, :] ./= Npairs
    end
    Fs .= sum(Fs_pp; dims=2)[:] / N
    return Fs, Fs_pp
end