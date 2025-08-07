"""
    find_overlap_function(s; a=0.5)

Calculates the overlap function for a simulation.

# Arguments
- `s`: The simulation.
- `a=0.5`: The cutoff distance.

# Returns
- `Fs`: The overlap function.
- `Fs_pp`: The overlap function per particle.
"""
function find_overlap_function(s; a=0.5)
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