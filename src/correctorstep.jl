function correctorstep!(A::Vector{Float64},Q::Vector{Float64},Aprev::Vector{Float64},
    Qprev::Vector{Float64},Fbarforward::Vector{Float64},Fbarbackward::Vector{Float64},
    Qforward::Vector{Float64},Aforward::Vector{Float64},Qbackward::Vector{Float64},Abackward::Vector{Float64},
    colsint::Vector{Int8},acolscor::Vector{Int8},qcolscor::Vector{Int8},
    h::Float64,k::Float64,diff::Float64,μ::Float64,rho::Float64)
    A[:] .= Aprev[colsint] .- h./k.*(Fbarforward[acolscor] .- Fbarbackward[acolscor]);
    Q[:] .= Qprev[colsint] .- h./k.*(Fbarforward[qcolscor] .- Fbarbackward[qcolscor]) .-
        h.*diff.*π.*μ./rho.*(0.5.*(Qforward./Aforward .+ Qbackward./Abackward));
end
