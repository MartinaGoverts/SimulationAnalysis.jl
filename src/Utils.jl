"""
    find_relaxation_time(t::AbstractVector, F::AbstractVector; threshold::Float64=exp(-1))

Calculates the relaxation time `τ` of a decaying function `F(t)`.

The relaxation time is defined as the time at which `F(t)` first drops below a certain `threshold`.
The function finds this time by performing linear interpolation between the two points straddling the threshold. The interpolation is linear in `F` and `log(t)`.

# Arguments
- `t::AbstractVector`: A vector of time points.
- `F::AbstractVector`: A vector of the function values `F(t)`, corresponding to the time points in `t`. This function should be monotonically decreasing.
- `threshold::Float64=exp(-1)`: The threshold value for `F(t)`. The relaxation time is `t` such that `F(t) = threshold`.

# Returns
- `τ::Float64`: The calculated relaxation time. Returns `0.0` if the function is already below the threshold at the first time point. Returns `Inf` if the function never drops below the threshold.
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