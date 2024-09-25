
"""
    find_correlation_function!(C, A, B, kspace, dt_array, t1_t2_pair_array)

Computes the correlation function between two lists of observables. The observables are assumed to be complex-valued. The correlation function is a complex number.
The correlation function is computed for a given k-space, a given set of time differences, and a given set of pairs of times. The correlation function is computed for each time difference and is stored in the array `C`.


# Arguments
- `C::Array{Complex{Float64}, 1}`: The array where the correlation function is stored.
- `A::Array{Complex{Float64}, 2}`: The first list of observables. The first index is the wavevector index and the second index is the time index.
- `B::Array{Complex{Float64}, 2}`: The second list of observables.
- `kspace::KSpace`: The k-space.
- `dt_array::Array{Float64, 1}`: The array of time differences.
- `t1_t2_pair_array::Array{Array{Int64, 2}, 1}`: The array of pairs of times. Each element of the array is a matrix with two columns. Each row of the matrix contains the indices of the times `t1` and `t2` for which the correlation function is computed.
- `kmin::Float64`: The minimum value of the wavevector.
- `kmax::Float64`: The maximum value of the wavevector.

# Returns
- `Nothing`: The correlation function is stored in the array `C`.

"""
function real_correlation_function!(C, ReA::Array{Float64, 2}, ImA::Array{Float64, 2}, ReB::Array{Float64, 2}, ImB::Array{Float64, 2}, 
                                    kspace::KSpace, dt_array, t1_t2_pair_array, kmin, kmax)
    # not normalized by N!
    klengths = kspace.k_lengths
    Nk_samples = 0
    _, Nk = size(ReA)
    @assert size(ReA) == size(ReB) == size(ImA) == size(ImB)
    @inbounds for ik = 1:Nk
        klength = klengths[ik]
        if !(kmin ≤ klength ≤ kmax)
            continue
        end
        Nk_samples += 1
        for iδt in eachindex(dt_array)
            pairs_idt = t1_t2_pair_array[iδt]
            Npairs = size(pairs_idt, 1)
            for ipair = 1:Npairs
                t1 = pairs_idt[ipair, 1]
                t2 = pairs_idt[ipair, 2]
                C[iδt] += ReB[t2, ik]*ReA[t1, ik] + ImB[t2, ik]*ImA[t1, ik]
            end
        end
    end
    for iδt in eachindex(dt_array)
        pairs_idt = t1_t2_pair_array[iδt]
        Npairs = size(pairs_idt, 1)
        C[iδt] /= Npairs*Nk_samples
    end
end


"""
    find_correlation_function!(AB, A, B, kspace::KSpace, dt_array, t1_t2_pair_array)

Computes the correlation function between two lists of observables. The observables are assumed to be complex-valued. The correlation function is a complex number.
The correlation function is averaged over the first dimension of A and B, a given set of time differences, and a given set of pairs of times. The correlation function is computed for each time difference and is stored in the array `AB`.

# Arguments
- `AB::Array{Complex{Float64}, 1}`: The array where the correlation function is stored.
- `A::Array{Complex{Float64}, 2}`: The first list of observables. The first index is the wavevector index and the second index is the time index.
- `B::Array{Complex{Float64}, 2}`: The second list of observables.
- `dt_array::Array{Float64, 1}`: The array of time differences.
- `t1_t2_pair_array::Array{Array{Int64, 2}, 1}`: The array of pairs of times. Each element of the array is a matrix with two columns. Each row of the matrix contains the indices of the times `t1` and `t2` for which the correlation function is computed.

# Returns
- `Nothing`: The correlation function is stored in the array `AB`.
"""
 function find_correlation_function!(AB, A, B, dt_array, t1_t2_pair_array)
    Nk = size(A, 1)
    @assert size(A) == size(B)
    for iδt in eachindex(dt_array)
         pairs_idt = t1_t2_pair_array[iδt]
         Npairs = size(pairs_idt, 1)
         for ipair = 1:Npairs
             t1 = pairs_idt[ipair, 1]
             t2 = pairs_idt[ipair, 2]
             AB_temp = 0.0+0.0im
             for ik = 1:Nk
                 AB_temp += conj(A[ik, t1])*B[ik, t2]
             end
             AB[iδt] += AB_temp
         end
     end
     for iδt in eachindex(dt_array)
         pairs_idt = t1_t2_pair_array[iδt]
         Npairs = size(pairs_idt, 1)
         AB[iδt] /= Nk*Npairs
     end
 end

"""
    find_correlation_function(A, B, s)

Computes the correlation function between two lists of observables. The observables are assumed to be complex-valued. The correlation function is a complex number.
The correlation function is averaged over the first dimension of A and B.

# Arguments
- `A::Array{Complex{Float64}, 2}`: The first list of observables.
- `B::Array{Complex{Float64}, 2}`: The second list of observables.
- `s::Simulation`: The simulation.

# Returns
- `C::Array{Complex{Float64}, 1}`: The correlation function.


"""
function find_correlation_function(A, B, s)
    dt_array = s.dt_array
    Ndt = length(dt_array)
    t1_t2_pair_array = s.t1_t2_pair_array
    C = zeros(ComplexF64, Ndt)
    find_correlation_function!(C, A, B, dt_array, t1_t2_pair_array)
    return C
end

"""
    find_correlation_matrix(observables_list1, observables_list2, kspace, s)

Computes the correlation matrix between two lists of observables. The observables are assumed to be complex-valued. The correlation matrix is a matrix of complex numbers.

# Arguments
- `observables_list1::Array{Array{Complex{Float64}, 2}, 1}`: The first list of observables.
- `observables_list2::Array{Array{Complex{Float64}, 2}, 1}`: The second list of observables.
- `s::Simulation`: The simulation.

# Returns
- `matrix::Array{Array{Complex{Float64}, 1}, 2}`: The correlation matrix.

"""
function find_correlation_matrix(observables_list1, observables_list2, s::Simulation)
    No1 = length(observables_list1)
    No2 = length(observables_list1)
    matrix = [zeros(ComplexF64, 0) for i=1:No1, j=1:No2] 
    for (i,observable_1) in enumerate(observables_list1)
        for (j,observable_2) in enumerate(observables_list2)
            matrix[i,j] = find_correlation_function(observable_1, observable_2, s)
        end
    end
    return matrix
end




"""
    real_static_correlation_function(ReA::Array{Float64, 2}, ImA::Array{Float64, 2}, ReB::Array{Float64, 2}, ImB::Array{Float64, 2}, kspace::KSpace, kmin, kmax)

Computes the real-valued static correlation function between two lists of observables. The observables are assumed to be complex-valued. The correlation function is a real number.

# Arguments
- `ReA::Array{Float64, 2}`: The real part of the first list of observables.
- `ImA::Array{Float64, 2}`: The imaginary part of the first list of observables.
- `ReB::Array{Float64, 2}`: The real part of the second list of observables.
- `ImB::Array{Float64, 2}`: The imaginary part of the second list of observables.
- `kspace::KSpace`: The k-space.
- `kmin::Float64`: The minimum value of the wavevector.
- `kmax::Float64`: The maximum value of the wavevector.

# Returns
- `C::Float64`: The correlation function.
"""
function real_static_correlation_function(ReA::Array{Float64, 2}, ImA::Array{Float64, 2}, ReB::Array{Float64, 2}, ImB::Array{Float64, 2}, 
                                    kspace::KSpace, kmin, kmax)
    # not normalized by N
    klengths = kspace.k_lengths
    Nk_samples = 0
    Nt, Nk = size(ReA)
    @assert size(ReA) == size(ReB) == size(ImA) == size(ImB)
    C = 0.0
    @inbounds for ik = 1:Nk
        klength = klengths[ik]
        if !(kmin ≤ klength ≤ kmax)
            continue
        end
        Nk_samples += 1
        for iδt in 1:Nt
            C += ReB[iδt, ik]*ReA[iδt, ik] + ImB[iδt, ik]*ImA[iδt, ik]
        end
    end
    C /= Nk_samples*Nt
    return C
end