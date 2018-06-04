function model1dhemo!(system::CVSystem,n::Int64,ID::Int64)

    if ~isempty(system.branches.children[ID])
        system.branches.A[ID][system.solverparams.JL,n+2] = (system.hemo.Ph/
            system.branches.beta[ID][end]+sqrt(system.branches.A0[ID][end]))^2;
        system.branches.W2[ID] = system.branches.W1[ID]-
            8*((system.branches.A[ID][system.solverparams.JL,n+2]*(system.branches.beta[ID][end]/
            (2*system.solverparams.rho))^2)^0.25-system.branches.c0[ID][end]);
        system.branches.Q[ID][system.solverparams.JL,n+2] = 0.5*
            system.branches.A[ID][system.solverparams.JL,n+2]*(system.branches.W1[ID]+
            system.branches.W2[ID])
    else
        system.branches.A[ID][system.solverparams.JL,n+2] = (system.hemo.Ph/
            system.branches.beta[ID][end]+sqrt(system.branches.A0[ID][end]))^2;
        W2 = system.branches.W1end[ID]-
            8*((system.branches.A[ID][system.solverparams.JL,n+2]*(system.branches.beta[ID][end]/
            (2*system.solverparams.rho))^2)^0.25-system.branches.c0[ID][end]);
        system.branches.Q[ID][system.solverparams.JL,n+2] = 0.5*
            system.branches.A[ID][system.solverparams.JL,n+2]*(system.branches.W1[ID]+
            W2)
    end

    system.hemo.Vloss += system.branches.Q[ID][system.solverparams.JL,n+2]*
        system.solverparams.h;

    if ~isempty(system.branches.children[ID]) # bifurcation downstream
        for i = 1:length(system.branches.children[ID])
            cID = system.branches.children[ID][i];
            system.branches.A[cID][1,n+2] = (system.hemo.Ph/system.branches.beta[cID][end]+
                sqrt(system.branches.A0[cID][end]))^2;
            system.branches.W1[cID] = system.branches.W2[ID]+8*((system.branches.A[cID][1,n+2]*
                (system.branches.beta[cID][end]/(2*system.solverparams.rho))^2)^0.25-
                system.branches.c0[cID][end]);
            system.branches.Q[cID][1,n+2] = 0.5*system.branches.A[cID][1,n+2]*
                (system.branches.W1[cID]+system.branches.W2[cID]);
            system.hemo.Vloss += -system.branches.Q[cID][1,n+2]*system.solverparams.h;
        end
    else # terminal compartments downstream
        system.branches.term[ID].V[n+2,1] = system.branches.term[ID].V[n+1,1]; # keep V constant (no pressure gradient to force out accumulating backflow)
        system.branches.term[ID].P[n+2,1] = system.hemo.Ph;
        system.hemo.Vloss += -system.branches.term[ID].Q[1,n+1]*system.solverparams.h;
    end

end
