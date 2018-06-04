function predictorstep!(Aforward::Vector{Float64},Qforward::Vector{Float64},
    Abackward::Vector{Float64},Qbackward::Vector{Float64},A::Vector{Float64},
    Q::Vector{Float64},h::Float64,k::Float64,Fp::Vector{Float64},colsint::Vector{Int8},
    acolspre::Vector{Int8},qcolspre::Vector{Int8},diff::Float64,μ::Float64,rho::Float64)
    Aforward .= (0.5.*(A[colsint.+1] .+ A[colsint]) .-
        0.5.*h./k.*(Fp[acolspre.+1] .- Fp[acolspre]));
    Qforward .= (0.5.*(Q[colsint.+1] .+ Q[colsint]) .- 0.5.*h./k.*(Fp[qcolspre.+1] .- Fp[qcolspre]) .-
        0.5.*h.*diff.*π.*μ./rho.*(0.5.*(Q[colsint.+1]./A[colsint.+1] .+ Q[colsint]./A[colsint])));
    Abackward .= (0.5.*(A[colsint] .+ A[colsint.-1]) .-
        0.5.*h./k.*(Fp[acolspre] .- Fp[acolspre.-1]));
    Qbackward .= (0.5.*(Q[colsint] .+ Q[colsint.-1]) .- 0.5.*h./k.*(Fp[qcolspre] .- Fp[qcolspre.-1]) .-
        0.5.*h.*diff.*π.*μ./rho.*(0.5.*(Q[colsint]./A[colsint] .+ Q[colsint.-1]./A[colsint.-1])));
end
