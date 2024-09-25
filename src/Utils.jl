function find_relaxation_time(t, F; threshold=exp(-1))
    logt = log.(t)
    for i in 1:length(t)
        if F[i] < threshold
            
            # linear interpolation on log t
            if i == 1
                return 0.0
            end
            slope = (logt[i] - logt[i-1])/(F[i] - F[i-1])
            return exp(logt[i-1] + slope*(threshold - F[i-1]))
        end
    end
    return Inf
end