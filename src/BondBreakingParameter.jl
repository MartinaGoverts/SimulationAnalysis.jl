
function CB_microkernel(neigh1::Vector{Int}, neigh2::Vector{Int})
    N_neig1 = length(neigh1)

    if N_neig1 == 0 
        return 1000.0
    end
    N_neig_still2 = 0
    for item in neigh1
        if item in neigh2
            N_neig_still2 += 1
        end
    end

    return N_neig_still2/N_neig1

end

function find_CB_per_particle(neighbourlists1, neighbourlists2, particle_i, dt_array, t1_t2_pair_array)
    Ndt = length(dt_array)
    CB = zeros(Ndt)

    for iδt in eachindex(dt_array)
        pairs_idt = t1_t2_pair_array[iδt]
        Npairs = size(pairs_idt, 1)
        C_temp = 0.0
        pairs_computed = 0
        for ipair = 1:Npairs
            t1 = pairs_idt[ipair, 1]
            neighbourlist1 = neighbourlists1[t1][particle_i]
            t2 = pairs_idt[ipair, 2]
            neighbourlist2 = neighbourlists2[t2][particle_i]

            CBnew = CB_microkernel(neighbourlist1, neighbourlist2)
            if CBnew < 100
                C_temp += CBnew
                pairs_computed += 1
            else
                # @show particle_i, t1
            end
        end
        CB[iδt] = C_temp
        pairs_idt = t1_t2_pair_array[iδt]
        Npairs = size(pairs_idt, 1)
        if pairs_computed == 0
            @show particle_i, iδt, Npairs
        end
        CB[iδt] /= pairs_computed
    end

    return CB
end

function find_CB(s, neighbourlists1, neighbourlists2)
    dt_arr = s.dt_array
    t1_t2_pair_array = s.t1_t2_pair_array
    N_particles = s.N
    CB = zeros(length(dt_arr), N_particles)
    for particle_i = 1:N_particles
        CB[:, particle_i] .= find_CB_per_particle(neighbourlists1, neighbourlists2, particle_i, dt_arr, t1_t2_pair_array)
    end
    return CB
end
