
"""
    find_diagonal_four_point_k_vector_set(kspace::KSpace3D, k1_bounds, k2_bounds, cosθ_bounds)    

    Finds the set of k vectors that satisfy the conditions given for k1, k2 and cosθ12.
    These are specified in 2-Tuples.

"""
function find_diagonal_four_point_k_vector_set(kspace, k1_bounds, k2_bounds, cosθ_bounds)
    #Initialize set
    k_vector_set = Tuple{Int64, Int64}[]

    k_array = kspace.k_array
    k_lengths = kspace.k_lengths
    Nk = length(k_lengths)

    # check bounds
    @assert k1_bounds[1] < k1_bounds[2]
    @assert k2_bounds[1] < k2_bounds[2]
    @assert cosθ_bounds[1] < cosθ_bounds[2]
        
    # find all suitable vectors
    ik1_set = [ik for ik = 1:Nk if k1_bounds[1] < k_lengths[ik] < k1_bounds[2]]
    ik2_set = [ik for ik = 1:Nk if k2_bounds[1] < k_lengths[ik] < k2_bounds[2]]

    for ik1 in ik1_set
        k1 = k_lengths[ik1]
        if k1 == 0
            continue
        end
        k1x = k_array[1, ik1]
        k1y = k_array[2, ik1]
        k1z = k_array[3, ik1]
        for ik2 in ik2_set
            if k1 == 0
                continue
            end
            if ik1 == ik2
                continue
            end
            k2 = k_lengths[ik2]
            k2x = k_array[1, ik2]
            k2y = k_array[2, ik2]
            k2z = k_array[3, ik2]

            k1dotk2 = k1x*k2x + k1y*k2y + k1z*k2z
            cosθ = k1dotk2 / (k1*k2)

            if cosθ_bounds[1] < cosθ < cosθ_bounds[2]
                push!(k_vector_set, (ik1, ik2))
            end
        end
    end
    kvecset = zeros(Int64, length(k_vector_set), 2)
    for ik in 1:length(k_vector_set)
        kvecset[ik, 1] = k_vector_set[ik][1]
        kvecset[ik, 2] = k_vector_set[ik][2]
    end
    return kvecset
end


function find_F4_diagonal(s::Simulation, kspace, k1_bounds, k2_bounds, cosθ_bounds)
    dt_array = s.dt_array
    t1_t2_pair_array = s.t1_t2_pair_array
    num_dt = length(dt_array)
    kvecset = find_diagonal_four_point_k_vector_set(kspace, k1_bounds, k2_bounds, cosθ_bounds)
    F4_arr = zeros(num_dt)
    F2F2_arr = zeros(num_dt)
    println("Calculating diagonal F4(k1=$((k1_bounds[2]+k1_bounds[1])/2), k2=$((k2_bounds[2]+k2_bounds[1])/2), cosθ=$((cosθ_bounds[2]+cosθ_bounds[1])/2), t)")

    @threads for idt in 1:num_dt
        pairs_idt = t1_t2_pair_array[idt]
        F4, F2F2 = find_F4_diagonal(kspace, kvecset, pairs_idt)
        F4_arr[idt] = F4
        F2F2_arr[idt] = F2F2
    end

    return F4_arr, F2F2_arr
end

function find_F4_diagonal(kspace, kvecset, pairs_idt)
    Reρkt = kspace.Reρkt
    Imρkt = kspace.Imρkt
    N = kspace.N
    F4 = 0.0
    F2F2 = 0.0
    Npairs = size(pairs_idt,1)
    for ipair = 1:Npairs
        t1 = pairs_idt[ipair,1]
        t2 = pairs_idt[ipair,2]
        F4new, F2F2new = _find_F4_diagonal(Reρkt, Imρkt, kvecset, N, t1, t2)
        F4 += F4new
        F2F2 += F2F2new
    end
    F4 /= Npairs
    F2F2 /= Npairs
    return F4, F2F2
end

function _find_F4_diagonal(Reρkt, Imρkt, kvecset, N, it1, it2)
    F4 = 0.0
    F2k1 = 0.0
    F2k2 = 0.0
    @inbounds for ivecs in 1:size(kvecset, 1)
        ik1 = kvecset[ivecs,1]
        ik2 = kvecset[ivecs,2]

        dF4, dF2k1, dF2k2 = _F4_kernel(Reρkt, Imρkt, ik1, ik2, it1, it2)

        F4 += dF4
        F2k1 += dF2k1
        F2k2 += dF2k2
    end
    F4 /= N^2*size(kvecset,1)
    F2k1 /= N*size(kvecset,1)
    F2k2 /= N*size(kvecset,1)
    return F4, F2k1*F2k2
end

Base.Base.@propagate_inbounds function _F4_kernel(Reρkt, Imρkt, ik1, ik2, it1, it2)
    Reρk1_0 = Reρkt[it1, ik1]
    Imρk1_0 = Imρkt[it1, ik1]
    Reρk2_0 = Reρkt[it1, ik2]
    Imρk2_0 = Imρkt[it1, ik2]
    Reρk1_t = Reρkt[it2, ik1]
    Imρk1_t = Imρkt[it2, ik1]
    Reρk2_t = Reρkt[it2, ik2]
    Imρk2_t = Imρkt[it2, ik2]

    ρk1ρk2_0_conjugated_real = Reρk1_0*Reρk2_0 - Imρk1_0*Imρk2_0
    ρk1ρk2_0_conjugated_imag = -Reρk1_0*Imρk2_0 - Reρk2_0*Imρk1_0

    ρk1ρk2_t_real = Reρk1_t*Reρk2_t - Imρk1_t*Imρk2_t
    ρk1ρk2_t_imag = Reρk1_t*Imρk2_t + Reρk2_t*Imρk1_t
    F4 = ρk1ρk2_0_conjugated_real*ρk1ρk2_t_real - ρk1ρk2_0_conjugated_imag*ρk1ρk2_t_imag
    F2k1 = Reρk1_0*Reρk1_t + Imρk1_0 * Imρk1_t
    F2k2 = Reρk2_0*Reρk2_t + Imρk2_0 * Imρk2_t
    return F4, F2k1, F2k2
end


function find_F4_diagonal_all_k(s, kspace, idt, Nksample, Ncostheta)
    k_lengths = kspace.k_lengths
    k_array = kspace.k_array
    dt_array = s.dt_array
    println("Finding F4(k1, k2, theta, t=$(dt_array[idt]))")

    t1_t2_pair_array = s.t1_t2_pair_array
    pairs_idt = t1_t2_pair_array[idt]
    Npairs = size(pairs_idt, 1)
    N = kspace.N
    Reρkt = kspace.Reρkt
    Imρkt = kspace.Imρkt
    Nk = length(k_lengths)
    kmax = maximum(k_lengths)/sqrt(3.0)
    dk = kmax/Nksample
    k_sample_array = collect(1:Nksample)*dk .- dk/2
    dcostheta = 2.0 / Ncostheta
    costheta_sample_array = collect(1:Ncostheta)*dcostheta .- dcostheta/2 .- 1.0
    F4_arr = zeros(Nksample, Nksample, Ncostheta)
    F2k1_arr = zeros(Nksample, Nksample, Ncostheta)
    F2k2_arr = zeros(Nksample, Nksample, Ncostheta)
    count = zeros(Int64, Nksample, Nksample, Ncostheta)

    @inbounds for ik1 in 1:Nk
        if ik1 % 10000 == 0
            println(ik1, "/", Nk)
        end
        k1 = k_lengths[ik1]
        if k1 == 0
            continue
        end
        k1x = k_array[1, ik1]
        k1y = k_array[2, ik1]
        k1z = k_array[3, ik1]
        i_sampled_k1 = floor(Int64, k1/dk) + 1
        if i_sampled_k1 > Nksample
            continue
        end
        for ik2 in 1:Nk
            if ik1 == ik2
                continue
            end

            k2 = k_lengths[ik2]
            if k2 == 0
                continue
            end
            i_sampled_k2 = floor(Int64, k2/dk) + 1
            if i_sampled_k2 > Nksample
                continue
            end
            k2x = k_array[1, ik2]
            k2y = k_array[2, ik2]
            k2z = k_array[3, ik2]
            if k2x == - k1x && k2y == -k1y && k2z == -k1z
                continue
            end
            k1dotk2 = k1x*k2x + k1y*k2y + k1z*k2z
            cosθ = k1dotk2 / (k1*k2)
            i_sampled_cosθ = floor(Int64, (cosθ + 1.0)/dcostheta) + 1
            if 0.999999 <= cosθ <= 1.000001
                i_sampled_cosθ = Ncostheta
            end
            if -0.999999 >= cosθ >= -1.000001
                i_sampled_cosθ = 1
            end
            dF4, dF2k1, dF2k2 = 0.0, 0.0, 0.0

            for ipair = 1:Npairs
                t1 = pairs_idt[ipair,1]
                t2 = pairs_idt[ipair,2]

                ddF4, ddF2k1, ddF2k2 = _F4_kernel(Reρkt, Imρkt, ik1, ik2, t1, t2)
                
                dF4 += ddF4
                dF2k1 += ddF2k1
                dF2k2 += ddF2k2
            end

            F4_arr[i_sampled_k1, i_sampled_k2, i_sampled_cosθ] += dF4 / Npairs
            F2k1_arr[i_sampled_k1, i_sampled_k2, i_sampled_cosθ] += dF2k1 / Npairs
            F2k2_arr[i_sampled_k1, i_sampled_k2, i_sampled_cosθ] += dF2k2 / Npairs
            count[i_sampled_k1, i_sampled_k2, i_sampled_cosθ] += 1
        end
    end
    F4_arr[count .> 0] ./= count[count .> 0]
    F2k1_arr[count .> 0] ./= count[count .> 0]
    F2k2_arr[count .> 0] ./= count[count .> 0]
    println("Evaluated F4 at $(sum(count)) sets of wave vectors")

    F4_arr ./= N^2
    F2k1_arr ./= N
    F2k2_arr ./= N
    F2F2_arr = F2k1_arr .* F2k2_arr
    return k_sample_array, costheta_sample_array, F4_arr, F2F2_arr
end


function find_F4_super_diagonal(s, kspace; dk=0.1, kmax=0.0)
    Nk = kspace.Nk
    k_lengths = kspace.k_lengths
    Reρkt = kspace.Reρkt
    Imρkt = kspace.Imρkt
    if kmax == 0.0
        kmax = maximum(k_lengths)/sqrt(2)
    end
    dt_array = s.dt_array
    t1_t2_pair_array = s.t1_t2_pair_array
    Ndt = length(dt_array)
    F4 = zeros(Nk, Ndt)
    F2 = zeros(Nk, Ndt)

    println("Calculating superdiagonal F4(k,t)")
    @time @threads for ik = 1:Nk
        if k_lengths[ik] > kmax + dk
            continue
        end

        for idt in eachindex(dt_array)
            pairs_idt = t1_t2_pair_array[idt]
            Npairs = size(pairs_idt,1)

            F4temp = 0.0
            F2temp = 0.0
            @turbo for ipair = 1:Npairs
                t1 = pairs_idt[ipair,1]
                t2 = pairs_idt[ipair,2]

                Reρ_0 = Reρkt[t1, ik]
                Imρ_0 = Imρkt[t1, ik]
                Reρ_t = Reρkt[t2, ik]
                Imρ_t = Imρkt[t2, ik]
        
                ρkρkt_real = Reρ_0*Reρ_t + Imρ_0*Imρ_t
                ρkρkt_imag = Reρ_0*Imρ_t - Imρ_0*Reρ_t

                F4temp += ρkρkt_real^2 - ρkρkt_imag^2
                F2temp += ρkρkt_real
            end
            F4[ik, idt] = F4temp
            F2[ik, idt] = F2temp
        end
    end
    for idt in eachindex(dt_array)
        pairs_idt = t1_t2_pair_array[idt]
        Npairs = size(pairs_idt,1)
        F4[:, idt] ./= Npairs*kspace.N^2
        F2[:, idt] ./= Npairs*kspace.N
    end

    k_sample_array, F4_binned, F2_binned = bin_F4_and_F2(F4, F2, kspace; dk=dk, kmax=kmax)
    return k_sample_array, F4_binned, F2_binned.^2
end

function bin_F4_and_F2(F4, F2, kspace; dk=0.1, kmax=0.0)
    Ndt = size(F2, 2)
    Nk = kspace.Nk
    k_lengths = kspace.k_lengths
    if kmax == 0.0
        kmax = maximum(k_lengths)/sqrt(2)
    end
    kmin = minimum(k_lengths)
    krange = kmax - kmin
    Nbins = floor(Int64, krange / dk)
    k_sample_array = collect(1:Nbins)*dk .- dk/2 .+ kmin
    F4_binned = zeros(Nbins, Ndt)
    F2_binned = zeros(Nbins, Ndt)
    counts = zeros(Nbins)
    println("Binning F4(k,t) into $Nbins bins with  width dk = $dk")

    @time for ik ∈ 1:Nk
        k_i = k_lengths[ik]
        for ik_sample ∈ 1:Nbins
            k_sample = k_sample_array[ik_sample]
            if abs(k_i - k_sample) < dk/2
                for it = 1:Ndt
                    F4_binned[ik_sample, it] += F4[ik, it]
                    F2_binned[ik_sample, it] += F2[ik, it]
                end
                counts[ik_sample] += 1
            end
        end
    end
    F4_binned[counts .> 0, :] ./= counts[counts .> 0]
    F2_binned[counts .> 0, :] ./= counts[counts .> 0]
    return k_sample_array, F4_binned, F2_binned
end