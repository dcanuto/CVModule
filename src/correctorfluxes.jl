function correctorfluxes!(Fbarforward::Vector{Float64},Fbarbackward::Vector{Float64},
    Qforward::Vector{Float64},Aforward::Vector{Float64},Qbackward::Vector{Float64},
    Abackward::Vector{Float64},beta::Float64,rho::Float64,
    acolscor::Vector{Int8},qcolscor::Vector{Int8})
    Fbarforward[acolscor] .= Qforward;
    Fbarforward[qcolscor] .= Qforward.^2./Aforward .+ (1./3).*beta./rho.*
        Aforward.^1.5;
        # system.branches.Fbarforward[i][system.solverparams.acolscor] .= (
        #     system.branches.Qforward[i]);
        # system.branches.Fbarforward[i][system.solverparams.qcolscor] .= (
        #     system.branches.Qforward[i].^2./system.branches.Aforward[i] .+
        #     (1./3).*system.branches.beta[i][end]./system.solverparams.rho.*
        #     system.branches.Aforward[i].^1.5);
    Fbarbackward[acolscor] .= Qbackward;
    Fbarbackward[qcolscor] .= Qbackward.^2./Abackward .+ (1./3).*beta./rho.*
        Abackward.^1.5;
        # system.branches.Fbarbackward[i][system.solverparams.acolscor] .= (
        #     system.branches.Qbackward[i]);
        # system.branches.Fbarbackward[i][system.solverparams.qcolscor] .= (
        #     system.branches.Qbackward[i].^2./system.branches.Abackward[i] .+
        #     (1./3).*system.branches.beta[i][end]./system.solverparams.rho.*
        #     system.branches.Abackward[i].^1.5);
end
