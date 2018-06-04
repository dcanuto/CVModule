function splitinvariants!(system::CVSystem,n::Int64)
    lf = Float64;
    lb = Float64;
    W1 = Vector{Float64}[];
    dW1 = Float64;
    W2 = Vector{Float64}[];
    dW2 = Float64;
    ID = Int64;

    for i = 1:length(system.branches.ID)
        if ~isempty(system.branches.children[i])
            # forward wave speed for parent's distal end
            lf = (system.branches.Q[i][system.solverparams.JL,n+1]/
                system.branches.A[i][system.solverparams.JL,n+1]+
                sqrt(0.5*system.branches.beta[i][end]/system.solverparams.rho)*
                system.branches.A[i][system.solverparams.JL,n+1]^0.25);
            # parent's right running invariant at current time
            W1 = ([system.branches.Q[i][system.solverparams.JL-1,n+1]/
                system.branches.A[i][system.solverparams.JL-1,n+1]+
                4*(sqrt(0.5*system.branches.beta[i][end]/
                system.solverparams.rho)*
                system.branches.A[i][system.solverparams.JL-1,n+1]^0.25 -
                system.branches.c0[i][end])]);
            push!(W1,system.branches.Q[i][system.solverparams.JL,n+1]/
                system.branches.A[i][system.solverparams.JL,n+1]+
                4*(sqrt(0.5*system.branches.beta[i][end]/
                system.solverparams.rho)*
                system.branches.A[i][system.solverparams.JL,n+1]^0.25 -
                system.branches.c0[i][end]));
            dW1 = (W1[2] - W1[1])/system.branches.k[i];
            dQ = (system.branches.Q[i][system.solverparams.JL,n+1]-
                system.branches.Q[i][system.solverparams.JL-1,n+1])/system.branches.k[i];
            dA = (system.branches.A[i][system.solverparams.JL,n+1]-
                system.branches.A[i][system.solverparams.JL-1,n+1])/system.branches.k[i];
            # advance parent's distal boundary W1 by linear interp.
            Qint = system.branches.Q[i][system.solverparams.JL,n+1] - dQ*lf*system.solverparams.h;
            Aint = system.branches.A[i][system.solverparams.JL,n+1] - dA*lf*system.solverparams.h;
            system.branches.W1[i] = W1[2] - dW1*lf*system.solverparams.h -
                system.solverparams.diffusioncoeff*pi*system.solverparams.mu/
                system.solverparams.rho*Qint/Aint^2*system.solverparams.h;
            # backward wave speed for each child's proximal end
            for j = 1:length(system.branches.children[i])
                ID = system.branches.children[i][j];
                lb = (system.branches.Q[ID][1,n+1]/system.branches.A[ID][1,n+1]-
                    sqrt(0.5*system.branches.beta[ID][end]/
                    system.solverparams.rho)*system.branches.A[ID][1,n+1]^0.25);
                # child's left-running invariant at current time
                W2 = ([system.branches.Q[ID][1,n+1]/
                    system.branches.A[ID][1,n+1]-
                    4*(sqrt(0.5*system.branches.beta[ID][end]/
                    system.solverparams.rho)*system.branches.A[ID][1,n+1]^0.25-
                    system.branches.c0[ID][end])]);
                push!(W2,system.branches.Q[ID][2,n+1]/
                    system.branches.A[ID][2,n+1]-
                    4*(sqrt(0.5*system.branches.beta[ID][end]/
                    system.solverparams.rho)*system.branches.A[ID][2,n+1]^0.25-
                    system.branches.c0[ID][end]))
                dW2 = (W2[2] - W2[1])/system.branches.k[ID];
                dQ = (system.branches.Q[ID][2,n+1]-
                    system.branches.Q[ID][1,n+1])/system.branches.k[ID];
                dA = (system.branches.A[ID][2,n+1]-
                    system.branches.A[ID][1,n+1])/system.branches.k[ID];
                # advance child's proximal boundary invariant
                Qint = system.branches.Q[ID][1,n+1] - dQ*lb*system.solverparams.h;
                Aint = system.branches.A[ID][1,n+1] - dA*lb*system.solverparams.h;
                system.branches.W2[ID] = W2[1] - dW2*lb*system.solverparams.h -
                    system.solverparams.diffusioncoeff*pi*system.solverparams.mu/
                    system.solverparams.rho*Qint/Aint^2*system.solverparams.h;;
            end
        end
    end
end
