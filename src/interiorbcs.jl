function interiorbcs!(system::CVSystem,n::Int64,split::Int64,hemoflag::String,Q::Vector{Float64},A::Vector{Float64})
    children = system.branches.children[split];
    if length(children) == 1
        f = CVModule.fsingle;
        J = CVModule.Jsingle;
    elseif length(children) == 2
        f = CVModule.fdouble;
        J = CVModule.Jdouble;
    elseif length(children) == 3
        f = CVModule.ftriple;
        J = CVModule.Jtriple;
    elseif length(children) == 4
        f = CVModule.fquad;
        J = CVModule.Jquad;
    end
    W = zeros(length(children)+1);
    beta = zeros(length(children)+1);
    A0 = zeros(length(children)+1);
    c0 = zeros(length(children)+1);
    # k = zeros(length(children)+1);
    # Q = zeros(length(children)+1);
    # Qo = zeros(length(children)+1);
    # A = zeros(length(children)+1);
    # Ao = zeros(length(children)+1);
    beta[1] = system.branches.beta[split][end];
    A0[1] = system.branches.A0[split][end];
    c0[1] = system.branches.c0[split][end];
    # k[1] = system.branches.k[split][end];
    # Q[1] = system.branches.Q[splits[i]][system.solverparams.JL,n+1];
    # Qo[1] = system.branches.Q[split][system.solverparams.JL-1,n+1];
    # A[1] = system.branches.A[splits[i]][system.solverparams.JL,n+1];
    # Ao[1] = system.branches.A[split][system.solverparams.JL-1,n+1];
    for j = 1:length(children)
        beta[j+1] = system.branches.beta[children[j]][end];
        A0[j+1] = system.branches.A0[children[j]][end];
        c0[j+1] = system.branches.c0[children[j]][end];
        # k[j+1] = system.branches.k[children[j]][end];
        # Q[j+1] = system.branches.Q[children[j]][1,n+1];
        # Qo[j+1] = system.branches.Q[children[j]][2,n+1];
        # A[j+1] = system.branches.A[children[j]][1,n+1];
        # Ao[j+1] = system.branches.A[children[j]][2,n+1];
    end
    # CVModule.splitinvariants!(W,Q,A,Qo,Ao,beta,c0,system.solverparams.rho,
    #     k,system.solverparams.h,system.solverparams.diffusioncoeff,
    #     system.solverparams.mu,children);
    # system.branches.W1[split] = W[1];
    # system.branches.W2[children] .= W[2:end];
    W[1] = system.branches.W1[split];
    W[2:end] .= system.branches.W2[children];
    Qnew = zeros(length(children)+1);
    Anew = zeros(length(children)+1);
    iters = Int64[0];
    if hemoflag == "no"
        CVModule.solvesplits!(iters,Qnew,Anew,children,Q,A,W,beta,A0,c0,
            system.solverparams.rho,f,J,system.solverparams.maxiter,
            system.solverparams.maxval,system.solverparams.epsJ,
            system.solverparams.epsN)
    elseif hemoflag == "yes"
        curloss = [0.];
        CVModule.solvesplits!(iters,Qnew,Anew,children,Q,A,W,beta,A0,c0,
            system.solverparams.rho,f,J,system.solverparams.maxiter,
            system.solverparams.maxval,system.solverparams.epsJ,
            system.solverparams.epsN,hemoflag,system.hemo.injured[split],
            system.solverparams.h,system.hemo.Ph,curloss);
        system.hemo.Vloss += curloss[1];
        if system.hemo.injured[split] == true
            system.branches.W2[split] = W[1];
            system.branches.W1[children] .= W[2:end];
        end
        CVModule.applytourniquet!(system,n); # turn off to allow continual bleeding
    end
    system.solverparams.totaliter += iters[1];
    system.branches.Q[split][system.solverparams.JL,n+2] = Qnew[1];
    system.branches.A[split][system.solverparams.JL,n+2] = Anew[1];
    for j = 1:length(children)
        system.branches.Q[children[j]][1,n+2] = Qnew[j+1];
        system.branches.A[children[j]][1,n+2] = Anew[j+1];
    end
end
