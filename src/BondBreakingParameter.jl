
"""
    CB_microkernel(neigh1::Vector{Int}, neigh2::Vector{Int}) -> Float64

Calculates the fraction of neighbors in `neigh1` that are also present in `neigh2`.

Returns 1000.0 if the initial neighbor list `neigh1` is empty.
"""
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

"""
    find_CB_per_particle(neighbourlists1, neighbourlists2, particle_i, dt_array, t1_t2_pair_array)

Calculates the bond-breaking correlation function for a single particle, averaged over time origins.
"""
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

"""
    find_CB(s::Simulation, neighbourlists1::Vector{<:Vector}, neighbourlists2::Vector{<:Vector})

Calculates the bond-breaking correlation function `C_B(t)` for all particles.

This function measures the fraction of neighbors that a particle at time `t_0` still has at a later time `t_0 + t`.
It is averaged over all particles and time origins `t_0`.

# Arguments
- `s::Simulation`: The simulation data, used for time information.
- `neighbourlists1::Vector{<:Vector}`: A vector of neighbor lists at the initial times (`t_0`).
- `neighbourlists2::Vector{<:Vector}`: A vector of neighbor lists at the final times (`t_0 + t`). Often, this is the same as `neighbourlists1`.

# Returns
- `CB::Matrix{Float64}`: A `(Ndt, N)` matrix where `CB[i, j]` is the bond-breaking correlation for particle `j` at time delay `dt_array[i]`.
"""
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


"""
    compute_smoothed_gaussian(r_array, i, j, r, σ, box_size, ::Val{3})

Computes a smoothed Gaussian function in 3D.
"""
function compute_smoothed_gaussian(r_array, i, j, r, σ, box_size, ::Val{3}) 
    ri = SVector{3, Float64}(r_array[1, i], r_array[2, i], r_array[3, i])
    rj = SVector{3, Float64}(r_array[1, j], r_array[2, j], r_array[3, j])
    V = box_size^3
    prefactor = π^(-3/2)*σ^3 / (V*r)
    gaussian = 0.0
    for n = (
        (-1, -1, -1),
        (-1, -1, 0),
        (-1, -1, 1),
        (-1, 0, -1),
        (-1, 0, 0), 
        (-1, 0, 1), 
        (-1, 1, -1),
        (-1, 1, 0), 
        (-1, 1, 1), 
        (0, -1, -1),
        (0, -1, 0), 
        (0, -1, 1), 
        (0, 0, -1), 
        (0, 0, 0),  
        (0, 0, 1),  
        (0, 1, -1), 
        (0, 1, 0),  
        (0, 1, 1),  
        (1, -1, -1),
        (1, -1, 0), 
        (1, -1, 1), 
        (1, 0, -1), 
        (1, 0, 0),  
        (1, 0, 1),  
        (1, 1, -1),
        (1, 1, 0),
        (1, 1, 1),
        )
        imagex, imagey, imagez = n
        r_ij = rj - ri + SVector{3, Float64}(imagex*box_size, imagey*box_size, imagez*box_size)
        r_ij_length = norm(r_ij, 2)
        if i == j && imagex == 0 && imagey == 0 && imagez == 0 # same particle, same image
            gaussian += 0.0 # 2*exp(-r^2/(4σ^2))/π   
        else 
            term1 = exp(-(r_ij_length - r)^2/(4σ^2))
            term2 = exp(-(r_ij_length + r)^2/(4σ^2))
            gaussian += prefactor*(term1 - term2)/r_ij_length
        end
    end
    return gaussian
end


"""
    compute_smoothed_gaussian(r_array, i, j, r, σ, box_size, ::Val{2})

Computes a smoothed Gaussian function in 2D.
"""
function compute_smoothed_gaussian(r_array, i, j, r, σ, box_size, ::Val{2})
    ri = SVector{2, Float64}(r_array[1, i], r_array[2, i])
    rj = SVector{2, Float64}(r_array[1, j], r_array[2, j])
    V = box_size^2
    prefactor = σ*π^(-1/2)/(2V)
    gaussian = 0.0
    for n = (
        (-1, -1),
        (-1, 0),
        (-1, 1),
        (0, -1),
        (0, 0),
        (0, 1),
        (1, -1),
        (1, 0),
        (1, 1),
        )
        imagex, imagey = n
        r_ij = rj - ri + SVector{2, Float64}(imagex*box_size, imagey*box_size)
        r_ij_length = norm(r_ij, 2)
        if i == j && imagex == 0 && imagey == 0 # same particle, same image
            gaussian += 0.0 # 2*exp(-r^2/(4σ^2))/π 
        else 
            term1 = exp(-(r_ij_length^2 + r^2)/(4σ^2))
            term2 = besseli0(r*r_ij_length/(2σ^2))
            gaussian += prefactor*(term1*term2)
        end
    end
    return gaussian
end




"""
    find_χBB_smoothed(s, neighbourlists1, neighbourlists2, r, σ, CB_mean)

Computes the smoothed bond-breaking susceptibility.
"""
function find_χBB_smoothed(s, neighbourlists1, neighbourlists2, r, σ, CB_mean)
    dims = size(s.r_array, 1)
    find_χBB_smoothed(s, neighbourlists1, neighbourlists2, r, σ, CB_mean, Val(dims))
end

function find_χBB_smoothed(s, neighbourlists1, neighbourlists2, r,  σ, CB_mean, valdims::Val{dims}) where dims
    dt_arr = s.dt_array
    t1_t2_pair_array = s.t1_t2_pair_array
    N_particles = s.N
    χBB = zeros(length(dt_arr))
    δCb = zeros(N_particles)
    r_array = s.r_array
    box_size = s.box_sizes[1]
    dt_array = s.dt_array

    for iδt in eachindex(dt_array)
        pairs_idt = t1_t2_pair_array[iδt]
        Npairs = 1#size(pairs_idt, 1)
        for ipair = 1:Npairs
            t1 = pairs_idt[ipair, 1]
            t2 = pairs_idt[ipair, 2]
            for particle_i = 1:N_particles
                neighbourlist1 = neighbourlists1[t1][particle_i]
                neighbourlist2 = neighbourlists2[t2][particle_i]
                Cbi = CB_microkernel(neighbourlist1, neighbourlist2)
                if Cbi < 10
                    δCb[particle_i] = Cbi - CB_mean[iδt]
                else
                    δCb[particle_i] = 0.0
                end
            end
            r_array_t2 = @views r_array[:, :, t2]
            for i in 1:N_particles
                for j in 1:N_particles
                    gaussian = compute_smoothed_gaussian(r_array_t2, i, j, r, σ, box_size, valdims)
                    χBB[iδt] += δCb[i]*δCb[j]*gaussian
                end
            end
        end
        χBB[iδt] /= Npairs
        
    end

    return χBB
end


"""
    find_chi_BB(s, neighbourlists1, neighbourlists2, r_bin_edges::AbstractRange, cb; verbose=true)

Computes the bond-breaking susceptibility.
"""
function find_chi_BB(s, neighbourlists1, neighbourlists2, r_bin_edges::AbstractRange, cb; verbose=true)
    Ndims = size(s.r_array, 1)
    if Ndims == 2
        return find_chi_BB_2D(s, neighbourlists1, neighbourlists2, r_bin_edges, cb, verbose)
    elseif Ndims == 3
        return find_chi_BB_3D(s, neighbourlists1, neighbourlists2, r_bin_edges, cb, verbose)
    else
        error("Only 2D and 3D are supported")
    end
end


"""
    find_chi_BB_3D(s, neighbourlists1, neighbourlists2, r_bin_edges::AbstractRange, cb, verbose)

Computes the bond-breaking susceptibility in 3D.
"""
function find_chi_BB_3D(s, neighbourlists1, neighbourlists2, r_bin_edges::AbstractRange, cb, verbose)
    box_size = s.box_sizes[1]
    r_bin_width = step(r_bin_edges)
    dt_array = s.dt_array
    t1_t2_pair_array = s.t1_t2_pair_array
    N_particles = s.N

    χBB = zeros(length(dt_array), length(r_bin_edges)-1)
    count_array = zeros(size(χBB)...)

    Cb_i_all = zeros(N_particles, length(dt_array))
    @threads for (iδt) in eachindex(dt_array)
        if verbose
            println("Computing χBB for iδt = $iδt")
        end
        cb_iδt = cb[iδt]
        pairs_idt = t1_t2_pair_array[iδt]
        Npairs = size(pairs_idt, 1)
        χB_temp = 0.0
        for ipair = 1:Npairs
            t1 = pairs_idt[ipair, 1]
            t2 = pairs_idt[ipair, 2]

            for particle_i = 1:N_particles
                neighbourlist1_i = neighbourlists1[t1][particle_i]
                neighbourlist2_i = neighbourlists2[t2][particle_i]
                Cb_i_all[particle_i, iδt] = CB_microkernel(neighbourlist1_i, neighbourlist2_i)
            end

            for particle_i = 1:N_particles
                xi  = s.r_array[1, particle_i, t2]
                yi  = s.r_array[2, particle_i, t2]
                zi  = s.r_array[3, particle_i, t2]
                CBnew_i = Cb_i_all[particle_i, iδt]
                
                for particle_j = particle_i+1:N_particles
                    xj  = s.r_array[1, particle_j, t2]
                    yj  = s.r_array[2, particle_j, t2]
                    zj  = s.r_array[3, particle_j, t2]
                    CBnew_j = Cb_i_all[particle_j, iδt]
                    if CBnew_i < 100 && CBnew_j < 100
                        χB_temp = (CBnew_i - cb_iδt)*(CBnew_j - cb_iδt)
                        dx = xi - xj
                        dx = dx - box_size*round(dx/box_size)
                        dy = yi - yj
                        dy = dy - box_size*round(dy/box_size)
                        dz = zi - zj
                        dz = dz - box_size*round(dz/box_size)
                        r = sqrt(dx*dx + dy*dy + dz*dz)
                        ir = ceil(Int, r/r_bin_width)
                        if 1 <= ir <= length(r_bin_edges)-1
                            count_array[iδt, ir] += 1
                            χBB[iδt, ir] += χB_temp
                        end
                    end
                end
            end
        end
    end

    r_bin_centers = (r_bin_edges[1:end-1] + r_bin_edges[2:end])/2
    return dt_array, r_bin_centers, count_array, χBB
end

