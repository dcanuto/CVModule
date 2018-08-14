function arteryodes!(A::Vector{Float64},Q::Vector{Float64},dU::Vector{Float64},
    Fp::Vector{Float64},colsint::Vector{Int8},acols::Vector{Int8},qcols::Vector{Int8},
    diff::Float64,μ::Float64,rho::Float64,ts::Float64,As::Float64,Qs::Float64,
    Binv::Matrix{Float64},C::Matrix{Float64})
    # solve system of derivatives
    dA = Binv*C*Fp[acols];
    dQ = Binv*C*Fp[qcols];
    # output only interior points
    dU[acols[1:end.-2]] .= -1.*dA[2:10];
    dU[qcols[1:end.-2].-2] .= -1.*dQ[2:10] .- diff.*π.*μ.*ts./(rho.*As).*Q[colsint]./A[colsint];
    # dU[acols[1:end.-2]] .= .-1./(xs.*k).*(Fp[acols[3:end]] .- Fp[acols[2:end.-1]]);
    # dU[qcols[1:end.-2].-2] .= .-1./(xs.*k).*(Fp[qcols[3:end]] .- Fp[qcols[2:end.-1]]) .-
    #     diff.*π.*μ.*ts./(rho.*As).*Q[colsint]./A[colsint];
end
