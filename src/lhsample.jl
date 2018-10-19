using MAT

"""
    lhsample(mins, maxs, n)
Randomly sample `n` vectors from the parallelogram defined
by `mins` and `maxs` using the Latin hypercube algorithm.
Returns `dims`×`n` matrix.
"""

function lhsample(mins::AbstractVector{T},
                  maxs::AbstractVector{T},
                  n::Integer) where T<:Number
    length(mins) == length(maxs) ||
        throw(DimensionMismatch("mins and maxs should have the same length"))
    all(xy -> xy[1] <= xy[2], zip(mins, maxs)) ||
        throw(ArgumentError("mins[i] should not exceed maxs[i]"))
    dims = length(mins)
    result = zeros(T, dims, n)
    cubedim = Vector{T}(n)
    @inbounds for i in 1:dims # sample ith parameter
        imin = mins[i]
        dimstep = (maxs[i] - imin) / n
        for j in 1:n # sample jth stratum
            cubedim[j] = imin + dimstep * (j - 1 + rand(T))
        end
        result[i, :] .= shuffle!(cubedim) # randomly permute samples
    end
    file = matopen("lhpoints.mat","w");
    write(file,"lhpoints",result)
    close(file)
    return result
end
