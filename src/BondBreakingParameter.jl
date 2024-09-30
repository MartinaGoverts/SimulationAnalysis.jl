
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


function find_points_unit_sphere(r, N_angles, ::Val{2})
    dtheta = 2π/N_angles
    theta_arr = range(dtheta/2, stop=2π-dtheta/2, length=N_angles)
    points = SVector{2, Float64}[]
    for i in 1:N_angles
        theta = theta_arr[i]
        point = SVector{2, Float64}(r*cos(theta), r*sin(theta))
        push!(points, point)
    end
    return points
end

function find_points_unit_sphere(r, N_angles, ::Val{3})
    dcostheta = 2/N_angles
    dphi = π/N_angles
    costheta_arr = range(-1+dcostheta/2, stop=1-dcostheta/2, length=N_angles)
    phi_arr = range(dphi/2, stop=π-dphi/2, length=N_angles)
    points = SVector{3, Float64}[]
    for costheta = costheta_arr
        theta = acos(costheta)
        for phi = phi_arr
            point = SVector{3, Float64}(r*cos(phi)*sin(theta), r*sin(phi)*sin(theta), r*cos(theta))
            push!(points, point)
        end
    end
    return points
end


function compute_smoothed_gaussian(r_array, i, j, points, σ, box_size)
    ri = SVector{3, Float64}(r_array[1, i], r_array[2, i], r_array[3, i])
    rj = SVector{3, Float64}(r_array[1, j], r_array[2, j], r_array[3, j])

    r_ij = rj - ri
    r_ij = r_ij - box_size*round.(r_ij./box_size)

    gaussian = 0.0
    for rvec in points
        x = rvec + r_ij
        x2 = dot(x, x)
        arg = -x2/(4σ^2)
        gaussian += exp(arg)
    end
    gaussian /= length(points)
    return gaussian
end


function find_χBB(s, neighbourlists1, neighbourlists2,neighbourlistsχBB, r, N_angles, σ, CB_mean)
    dims = size(s.r_array, 1)
    find_χBB(s, neighbourlists1, neighbourlists2, neighbourlistsχBB, r, N_angles, σ, CB_mean, Val(dims))
end

function find_χBB(s, neighbourlists1, neighbourlists2, neighbourlistsχBB, r, N_angles, σ, CB_mean, valdims::Val{dims}) where dims
    dt_arr = s.dt_array
    t1_t2_pair_array = s.t1_t2_pair_array
    N_particles = s.N
    χBB = zeros(length(dt_arr))
    δCb = zeros(N_particles)
    r_array = s.r_array
    box_size = s.box_sizes[1]
    dt_array = s.dt_array
    points = find_points_unit_sphere(r, N_angles, valdims)

    for iδt in eachindex(dt_array)
        @show iδt
        pairs_idt = t1_t2_pair_array[iδt]
        Npairs = size(pairs_idt, 1)
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
            neighlist = neighbourlistsχBB[t2]
            for i in 1:N_particles
                for j in neighlist[i]
                    gaussian = compute_smoothed_gaussian(r_array_t2, i, j, points, σ, box_size)
                    χBB[iδt] += δCb[i]*δCb[j]*gaussian
                end
            end
        end
        χBB[iδt] *=  σ*sqrt(π)/(2prod(s.box_sizes)*Npairs)
    end

    return χBB
end