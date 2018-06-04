function predictorfluxes!(Fp::Vector{Float64},Q::Vector{Float64},A::Vector{Float64},
    beta::Float64,rho::Float64,acols::Vector{Int8},qcols::Vector{Int8})
    Fp[acols] .= Q;
    Fp[qcols] .= Q.^2./A .+ (1./3).*beta./rho.*A.^1.5;
end
