function model1dhemo!(yout::Vector{Float64},children::Vector{Int64},Ph::Float64,
    beta::Vector{Float64},rho::Float64,c0::Vector{Float64},A0::Vector{Float64},
    h::Float64,W1::Float64,W2 = [0.],V = 0.,Q = 0.)

    if ~isempty(children) # bifurcation downstream
        yout[1] = (Ph./beta[1] .+ sqrt.(A0[1])).^2;
        yout[2] = W1 .- 8.*((yout[1].*(0.5.*beta[1]./rho).^2).^0.25 .- c0[1]);
        yout[3] = 0.5.*yout[1].*(W1 .+ yout[2]);
        # system.branches.A[ID][system.solverparams.JL,n+2] = (system.hemo.Ph/
        #     system.branches.beta[ID][end]+sqrt(system.branches.A0[ID][end]))^2;
        # system.branches.W2[ID] = system.branches.W1[ID]-
        #     8*((system.branches.A[ID][system.solverparams.JL,n+2]*(system.branches.beta[ID][end]/
        #     (2*system.solverparams.rho))^2)^0.25-system.branches.c0[ID][end]);
        # system.branches.Q[ID][system.solverparams.JL,n+2] = 0.5*
        #     system.branches.A[ID][system.solverparams.JL,n+2]*(system.branches.W1[ID]+
        #     system.branches.W2[ID])
        yout[4] = yout[3].*h;
    else # terminal downstream
        yout[1] = (Ph./beta[1] .+ sqrt.(A0[1])).^2;
        W2 = W1 .- 8.*((yout[1].*(0.5.*beta[1]./rho).^2).^0.25 .- c0[1]);
        yout[2] = 0.5.*yout[1].*(W1 .+ W2);
        # system.branches.A[ID][system.solverparams.JL,n+2] = (system.hemo.Ph/
        #     system.branches.beta[ID][end]+sqrt(system.branches.A0[ID][end]))^2;
        # W2 = system.branches.W1end[ID]-
        #     8*((system.branches.A[ID][system.solverparams.JL,n+2]*(system.branches.beta[ID][end]/
        #     (2*system.solverparams.rho))^2)^0.25-system.branches.c0[ID][end]);
        # system.branches.Q[ID][system.solverparams.JL,n+2] = 0.5*
        #     system.branches.A[ID][system.solverparams.JL,n+2]*(system.branches.W1[ID]+
        #     W2)
        yout[3] = yout[2].*h;
    end


    # system.hemo.Vloss += system.branches.Q[ID][system.solverparams.JL,n+2]*
    #     system.solverparams.h;

    if ~isempty(children) # bifurcation downstream
        for i = 1:length(children)
            cID = children[i];
            yout[5+3*(i-1)] = (Ph./beta[i+1] .+ sqrt.(A0[i+1])).^2;
            yout[6+3*(i-1)] = W2[i] + 8.*((yout[5+3*(i-1)].*(0.5.*beta[i+1]./rho).^2)^0.25 .- c0[i+1]);
            yout[7+3*(i-1)] = 0.5.*yout[5+3*(i-1)].*(yout[6+3*(i-1)] + W2[i]);
            yout[4] += yout[7+3*(i-1)].*h;
            # system.branches.A[cID][1,n+2] = (system.hemo.Ph/system.branches.beta[cID][end]+
            #     sqrt(system.branches.A0[cID][end]))^2;
            # system.branches.W1[cID] = system.branches.W2[cID]+8*((system.branches.A[cID][1,n+2]*
            #     (system.branches.beta[cID][end]/(2*system.solverparams.rho))^2)^0.25-
            #     system.branches.c0[cID][end]);
            # system.branches.Q[cID][1,n+2] = 0.5*system.branches.A[cID][1,n+2]*
            #     (system.branches.W1[cID]+system.branches.W2[cID]);
            # system.hemo.Vloss += -system.branches.Q[cID][1,n+2]*system.solverparams.h;
        end
    else # terminal compartments downstream
        yout[4] = V;
        yout[5] = Ph;
        yout[3] += -Q.*h;
        # system.branches.term[ID].V[n+2,1] = system.branches.term[ID].V[n+1,1]; # keep V constant (no pressure gradient to force out accumulating backflow)
        # system.branches.term[ID].P[n+2,1] = system.hemo.Ph;
        # system.hemo.Vloss += -system.branches.term[ID].Q[1,n+1]*system.solverparams.h;
    end

end
