
function find_relative_distance_neighborlists(s, rc; ζ = 0.2)
    r_array = s.r_array
    D_array = s.D_array
    Nt = size(r_array, 3)
    N = size(r_array, 2)
    rc2 = rc^2
    neighbourlists = Vector{Vector{Vector{Int}}}(undef, Nt)
    x = rand(size(r_array, 1), N)*s.box_sizes[1]
    box = CellListMap.Box(s.box_sizes, maximum(s.D_array)*3)
    celllist = CellList(x, box, parallel=false)
    auxilliary_struct = CellListMap.AuxThreaded(celllist)

    function push_pair!(i, j, d2, neighbourlist_t, D_array, rc2, ζ)
        i == j && return
        Di = D_array[i]
        Dj = D_array[j]
        mean_d = (Di + Dj)/2.0 * (1 - ζ * abs(Di - Dj))
        mean_d_squared = mean_d*mean_d 

        if d2 / mean_d_squared < rc2
            push!(neighbourlist_t[i], j)
            push!(neighbourlist_t[j], i)

        end
        return neighbourlist_t
    end

    for t = 1:Nt
        neighbourlist_t = Vector{Vector{Int}}(undef, N)
        for i = 1:N
            neighbourlist_t[i] = Vector{Int}()
            sizehint!(neighbourlist_t[i], 8)
        end
        celllist = @views UpdateCellList!(r_array[:, :, t], box, celllist, auxilliary_struct, parallel=false)

        map_pairwise!(
            (x,y,i,j,d2,neighbourlist_t) -> push_pair!(i, j, d2, neighbourlist_t, D_array, rc2, ζ),
            neighbourlist_t,
            box,
            celllist,
            parallel=false
        )

        neighbourlists[t] = neighbourlist_t
    end
    return neighbourlists
end


function find_absolute_distance_neighborlists(s, rc)
    r_array = s.r_array
    D_array = s.D_array
    Nt = size(r_array, 3)
    N = size(r_array, 2)
    rc2 = rc^2
    neighbourlists = Vector{Vector{Vector{Int}}}(undef, Nt)
    x = rand(size(r_array, 1), N)*s.box_sizes[1]
    box = CellListMap.Box(s.box_sizes, min(s.box_sizes[1], rc))
    celllist = CellList(x, box, parallel=false)
    auxilliary_struct = CellListMap.AuxThreaded(celllist)

    function push_pair!(i, j, d2, neighbourlist_t, rc2)
        i == j && return
        if d2  < rc2
            push!(neighbourlist_t[i], j)
            push!(neighbourlist_t[j], i)
        end
        return neighbourlist_t
    end

    for t = 1:Nt
        neighbourlist_t = Vector{Vector{Int}}(undef, N)
        for i = 1:N
            neighbourlist_t[i] = Vector{Int}()
            sizehint!(neighbourlist_t[i], 8)
        end
        celllist = @views UpdateCellList!(r_array[:, :, t], box, celllist, auxilliary_struct, parallel=false)

        map_pairwise!(
            (x,y,i,j,d2,neighbourlist_t) -> push_pair!(i, j, d2, neighbourlist_t, rc2),
            neighbourlist_t,
            box,
            celllist,
            parallel=false
        )

        neighbourlists[t] = neighbourlist_t
    end
    return neighbourlists
end



function fill_r_array_with_images!(r_array_with_images, is_image_of, r_array, it, box_size, max_distance_from_boundary)
    Ndims, N, Nt = size(r_array)
    @assert it <= Nt
    index = 0
    
    for i in [0, 1, -1]
        for j in [0, 1, -1]
            if Ndims == 2
                for particle in 1:N
                    newpos = r_array[:, particle, it] .+ [i*box_size, j*box_size]
                    condition1 = newpos[1] < -max_distance_from_boundary
                    condition2 = newpos[1] > box_size + max_distance_from_boundary
                    condition3 = newpos[2] < -max_distance_from_boundary
                    condition4 = newpos[2] > box_size + max_distance_from_boundary
                    if condition1 || condition2 || condition3 || condition4
                        continue
                    end

                    index += 1

                    r_array_with_images[index] = SVector{2, Float64}(newpos)
                    is_image_of[index] = particle
                end
            elseif Ndims == 3
                for k in [0, 1, -1]
                    for particle in 1:N

                        newpos = r_array[:, particle, it] .+ [i*box_size, j*box_size, k*box_size]

                        condition1 = newpos[1] < -max_distance_from_boundary
                        condition2 = newpos[1] > box_size + max_distance_from_boundary
                        condition3 = newpos[2] < -max_distance_from_boundary
                        condition4 = newpos[2] > box_size + max_distance_from_boundary
                        condition5 = newpos[3] < -max_distance_from_boundary
                        condition6 = newpos[3] > box_size + max_distance_from_boundary

                        if condition1 || condition2 || condition3 || condition4 || condition5 || condition6
                            continue
                        end
                        index += 1

                        r_array_with_images[index] = SVector{3, Float64}(newpos)
                        is_image_of[index] = particle
                    end
                end
            end
        end
    end
    return index
end


function new_delaunay_facets(dhull)
    hull = dhull.hull
    maxlift = maximum(last, hull.pts)

    return filter(hull.facets.arr) do facet
        D = length(first(hull.pts))
        above_pt = sum(hull.pts[i] for i in facet.plane.point_indices) / D
        above_pt = setindex(above_pt, above_pt[end] + 2maxlift, D)
        Quickhull.hyperplane_dist(facet.plane, above_pt, hull.pts) < 0
    end
end

function new_facets(dhull)
    return Quickhull.mappedarray(f -> Quickhull.NgonFace(f.plane.point_indices), new_delaunay_facets(dhull))
end




"""
    find_neighborlists_voronoi(s::Simulation; max_distance_from_boundary=3.0) where Ndims

Computes the Voronoi neighborlists for a simulation.

# Arguments
- `s::Simulation`: The simulation.
- `max_distance_from_boundary::Float64`: The maximum distance from the boundary to take periodic images into account.

# Returns
- `neighbourlists::Vector{Vector{Vector{Int}}}`: The Voronoi neighborlists.

"""
function find_voronoi_neighborlists(s; max_distance_from_boundary=3.0, verbose=true)
    dims = size(s.r_array, 1)
    if dims == 2
        return find_voronoi_neighborlists(s, Val{2}(), max_distance_from_boundary, verbose)
    elseif dims == 3
        return find_voronoi_neighborlists(s, Val{3}(), max_distance_from_boundary, verbose)
    else
        throw(ArgumentError("Only 2D and 3D simulations are supported."))
    end
end


function find_voronoi_neighborlists(s, ::Val{Ndims}, max_distance_from_boundary, verbose) where Ndims
    r_array = s.r_array
    Nt = size(r_array, 3)
    N = size(r_array, 2)
    @assert Ndims == 2 || Ndims == 3
    @assert size(r_array, 1) == Ndims
    # Ndims = size(r_array, 1)
    @assert Ndims == 2 || Ndims == 3

    neighbourlists = Vector{Vector{Vector{Int}}}(undef, Nt)
    box_size = s.box_sizes[1]
    r_array_with_images = Vector{SVector{Ndims, Float64}}(undef, 10N)
    is_image_of = zeros(Int, 10N)
    for t = 1:Nt
        if verbose && t % 100 == 0
            println("$t / $Nt")
        end
        num_particles_with_images = fill_r_array_with_images!(r_array_with_images, is_image_of, r_array, t, box_size, max_distance_from_boundary)

        tri = Quickhull.delaunay(@views(r_array_with_images[1:num_particles_with_images]))
        facs = new_facets(tri)

        neighbourlist_t = [Set{Int}() for _ in 1:N]

        for triangle in facs
            triangle_points = triangle.data
            particle1 = triangle_points[1]
            particle1_original = is_image_of[particle1]
            particle2 = triangle_points[2]
            particle2_original = is_image_of[particle2]
            particle3 = triangle_points[3]
            particle3_original = is_image_of[particle3]

    
            if particle1_original == particle1
                push!(neighbourlist_t[particle1], particle2_original)
                push!(neighbourlist_t[particle1], particle3_original)
            end
    
            if particle2_original == particle2
                push!(neighbourlist_t[particle2], particle1_original)
                push!(neighbourlist_t[particle2], particle3_original)
            end
    
            if particle3_original == particle3
                push!(neighbourlist_t[particle3], particle1_original)
                push!(neighbourlist_t[particle3], particle2_original)
            end
    

            if Ndims == 3
                particle4 = triangle_points[4]
                particle4_original = is_image_of[particle4]
                if particle1_original == particle1
                    push!(neighbourlist_t[particle1], particle4_original)
                end
                if particle2_original == particle2
                    push!(neighbourlist_t[particle2], particle4_original)
                end
                if particle3_original == particle3
                    push!(neighbourlist_t[particle3], particle4_original)
                end
            
                if particle4_original == particle4
                    push!(neighbourlist_t[particle4], particle1_original)
                    push!(neighbourlist_t[particle4], particle2_original)
                    push!(neighbourlist_t[particle4], particle3_original)
                end
            end
        end

        neighbourlists[t] = collect.(neighbourlist_t)
    end
    return neighbourlists
end
