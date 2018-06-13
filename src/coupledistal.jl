function coupledistal!(system::CVSystem,n::Int64,terms::Vector{Int64},hemoflag="no")
    # W1 at next time step
    for i = 1:length(terms)
        ret1 = Float64[0];
        CVModule.endinvariants!(ret1,system.branches.Q[terms[i]][system.solverparams.JL,n+1],
            system.branches.A[terms[i]][system.solverparams.JL,n+1],system.branches.Q[terms[i]][system.solverparams.JL-1,n+1],
            system.branches.A[terms[i]][system.solverparams.JL-1,n+1],system.branches.c0[terms[i]][end],
            system.branches.beta[terms[i]][end],system.solverparams.rho,system.solverparams.mu,
            system.solverparams.diffusioncoeff,system.branches.k[terms[i]],system.solverparams.h);
        system.branches.W1end[terms[i]] = ret1[1];
    end
    if hemoflag == "no"
        for i = 1:length(terms)
            yout = zeros(4);
            Vs = system.branches.term[terms[i]].V[n+1,1];
            vs = system.branches.c0[terms[i]][end];
            ts = system.branches.k[terms[i]]/vs;
            iters = Int64[0];
            CVModule.newtondist!(yout,iters,Vs,vs,ts,system.branches.term[terms[i]].V[n+1,1],
                system.branches.term[terms[i]].C[1],system.branches.term[terms[i]].Q[n+1,1],
                system.branches.term[terms[i]].V0[1],system.branches.beta[terms[i]][end],
                system.branches.c0[terms[i]][end],system.branches.A0[terms[i]][end],
                system.branches.W1end[terms[i]],system.solverparams.rho,CVModule.fdist,
                CVModule.Jdist,system.solverparams.maxiter,system.solverparams.epsJ,
                system.solverparams.epsN,system.solverparams.maxval,system.solverparams.h);
            system.solverparams.totaliter += iters[1];
            system.branches.term[terms[i]].V[n+2,1] = yout[1];
            system.branches.term[terms[i]].P[n+2,1] = yout[2];
            system.branches.A[terms[i]][system.solverparams.JL,n+2] = yout[3];
            system.branches.Q[terms[i]][system.solverparams.JL,n+2] = yout[4];
        end
    elseif hemoflag == "yes"
        for i = 1:length(terms)
            if system.hemo.injured[terms[i]] == false
                yout = zeros(4);
                Vs = system.branches.term[terms[i]].V[n+1,1];
                vs = system.branches.c0[terms[i]][end];
                ts = system.branches.k[terms[i]]/vs;
                iters = Int64[0];
                CVModule.newtondist!(yout,iters,Vs,vs,ts,system.branches.term[terms[i]].V[n+1,1],
                    system.branches.term[terms[i]].C[1],system.branches.term[terms[i]].Q[n+1,1],
                    system.branches.term[terms[i]].V0[1],system.branches.beta[terms[i]][end],
                    system.branches.c0[terms[i]][end],system.branches.A0[terms[i]][end],
                    system.branches.W1end[terms[i]],system.solverparams.rho,CVModule.fdist,
                    CVModule.Jdist,system.solverparams.maxiter,system.solverparams.epsJ,
                    system.solverparams.epsN,system.solverparams.maxval,system.solverparams.h);
                system.solverparams.totaliter += iters[1];
                system.branches.term[terms[i]].V[n+2,1] = yout[1];
                system.branches.term[terms[i]].P[n+2,1] = yout[2];
                system.branches.A[terms[i]][system.solverparams.JL,n+2] = yout[3];
                system.branches.Q[terms[i]][system.solverparams.JL,n+2] = yout[4];
            elseif system.hemo.injured[terms[i]] == true
                children = system.branches.children[terms[i]];
                if ~isempty(children) # bifurncation downstream
                    beta = zeros(length(children[terms[i]])+1);
                    A0 = zeros(length(children[terms[i]])+1);
                    c0 = zeros(length(children[terms[i]])+1);
                    W2 = zeros(length(children[terms[i]]));
                    yout = zeros(4+3*length(children));
                    beta[1] = system.branches.beta[terms[i]][end];
                    A0[1] = system.branches.A0[terms[i]][end];
                    c0[1] = system.branches.c0[terms[i]][end];
                    for i = 1:length(children)
                        beta[i] = system.branches.beta[children[i]][end];
                        A0[i] = system.branches.A0[children[i]][end];
                        c0[i] = system.branches.c0[children[i]][end];
                        W2[i] = system.branches.W2[children[i]];
                    end
                    CVModule.model1dhemo!(yout,children,system.hemo.Ph,beta,
                        system.solverparams.rho,c0,A0,system.solverparams.h,
                        system.branches.W1[terms[i]],W2);
                    system.branches.A[terms[i]][system.solverparams.JL,n+2] = yout[1];
                    system.branches.W2[terms[i]] = yout[2];
                    system.branches.Q[terms[i]][system.solverparams.JL,n+2] = yout[3];
                    system.hemo.Vloss += yout[4];
                    for i = 1:length(children)
                        system.branches.A[children[i]][1,n+2] = yout[5+3*(i-1)];
                        system.branches.W1[children[i]] = yout[6+3*(i-1)];
                        system.branches.Q[children[i]][1,n+2] = yout[7+3*(i-1)];
                    end
                else # terminal downstream
                    yout = zeros(5);
                    beta = [system.branches.beta[terms[i]][end]];
                    A0 = [system.branches.A0[terms[i]][end]];
                    c0 = [system.branches.c0[terms[i]][end]];
                    CVModule.model1dhemo!(yout,children,system.hemo.Ph,beta,
                        system.solverparams.rho,c0,A0,system.solverparams.h,
                        system.branches.W1end[terms[i]],0.,system.branches.term[terms[i]].V[n+1,1],
                        system.branches.term[terms[i]].Q[n+1,1]);
                    system.branches.A[terms[i]][system.solverparams.JL,n+2] = yout[1];
                    system.branches.Q[terms[i]][system.solverparams.JL,n+2] = yout[2];
                    system.hemo.Vloss += yout[3];
                    system.branches.term[terms[i]].V[n+2,1] = yout[4];
                    system.branches.term[terms[i]].P[n+2,1] = yout[5];
                end
            end
        end
    end
end
