"""
    find_relaxation_time(t, F; threshold=exp(-1))

Calculates the relaxation time of a function `F(t)`.

The relaxation time is defined as the time at which `F(t)` drops below a certain `threshold`.
The function uses linear interpolation on a logarithmic time scale to find the relaxation time.

# Arguments
- `t`: The time points.
- `F`: The function values.
- `threshold=exp(-1)`: The threshold.

# Returns
- The relaxation time.
"""
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