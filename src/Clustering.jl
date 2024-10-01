

function find_cluster!(cluster, particle_i, neighbourlists, mobile_particles)

    for neighbour in neighbourlists[particle_i]
        if neighbour in mobile_particles
            push!(cluster, neighbour)
            delete!(mobile_particles, neighbour)
            find_cluster!(cluster, neighbour, neighbourlists, mobile_particles)
        end
    end

    return cluster
end


"""
    find_mobile_clusters(mobile_particles, neighbourlists)

Find clusters of mobile particles in a system. A cluster is defined as a set of particles that are connected by a bond. The function returns a vector of sets, where each set contains the indices of the particles in a cluster.

# Arguments
- `mobile_particles::Set{Int}`: The indices of the mobile particles. 
- `neighbourlists::Vector{Vector{Vector{Int}}}`: The neighbourlists of the particles.

# Returns
- `clusters::Vector{Set{Int}}`: The clusters of mobile particles.

"""
function find_mobile_clusters(mobile_particles::Set{Int}, neighbourlists)
    clusters = Vector{Set{Int}}()
    while !isempty(mobile_particles)
        particle_i = first(mobile_particles)
        delete!(mobile_particles, particle_i)
        cluster = Set{Int}(particle_i)
        find_cluster!(cluster, particle_i, neighbourlists, mobile_particles)
        push!(clusters, cluster)
    end
    return clusters
end