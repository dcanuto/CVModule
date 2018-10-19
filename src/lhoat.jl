using MAT

function lhoat(lhpoints::AbstractMatrix{T},
               f::AbstractVector{T}) where T <: Number
    length(f) == length(lhpoints[:,1]) ||
        throw(DimensionMismatch["need as many fractional changes as parameters"])
    nparams = length(lhpoints[:,1])
    npoints = length(lhpoints[1,:])
    for i in 1:npoints
        result = zeros(T,nparams,nparams)
        sample = Vector{T}(nparams)
        sample = copy(lhpoints[:,i])
        for j in 1:nparams
            sample[j] *= (1+f[j])
            result[:,j] = copy(sample)
            sample[j] /= (1+f[j])
        end
        file = matopen("lhoat_$(i).mat","w")
        write(file,"lhoat_$(i)",result)
        close(file)
    end
end
