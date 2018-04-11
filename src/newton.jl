function newton!(system::CVSystem,n::Int64,parentID::Int64)
    children = system.branches.children[parentID];
    numchildren = length(children);

    # initial guess state vector (flow rate/area pairs)
    x0 = zeros(2*numchildren+2);
    x0[1] = system.branches.Q[parentID][n+1,system.solverparams.JL];
    x0[2] = system.branches.A[parentID][n+1,system.solverparams.JL];
    for i = 1:numchildren
        ID = system.branches.children[parentID][i];
        x0[2*i+1] = system.branches.Q[ID][n+1,1];
        x0[2*i+2] = system.branches.A[ID][n+1,1];
    end

    # Newton iteration system of eqs. Jacobian based on junction type
    if numchildren == 1
        f = CVModule.fsingle;
        J = CVModule.Jsingle;
    elseif numchildren == 2
        f = CVModule.fdouble;
        J = CVModule.Jdouble;
    elseif numchildren == 3
        f = CVModule.ftriple;
        J = CVModule.Jtriple;
    elseif numchildren == 4
        f = CVModule.fquad;
        J = CVModule.Jquad;
    end

    xx = x0;
    x = zeros(length(x0));
    N = 1; # iteration counter

    # Newton iterations
    while N <= system.solverparams.maxiter
        # determine Jacobian, check invertibility
        JJ = J(xx,system,parentID,children);
        if abs(det(JJ)) < system.solverparams.epsJ
            error("Newton Jacobian is singular.");
        end
        # update state vector
        xn = xx - inv(JJ)*f(xx,system,parentID,children);
        # check for negative area (non-physical values)
        if xn[2] < 0
            error("Parent split area is negative. Parent ID: $parentID. Time step: $n.")
        end
        for i = 1:numchildren
            if xn[2*i+2] < 0
                error("Child split area is negative. Child ID: $(children[i]).
                    Parent ID: $parentID. Time step: $n.")
            end
        end
        # check if sufficiently close to root
        if norm(f(xn,system,parentID,children)) <= system.solverparams.epsN
            x = xn;
            # println(x)
            system.solverparams.totaliter+=N;
            break
        end
        # check for divergence
        if norm(f(xn,system,parentID,children)) >= system.solverparams.maxval
            system.solverparams.totaliter+=N;
            error("Newton iteration diverged.");
        end
        N+=1;
        xx = xn;
        if N == system.solverparams.maxiter
            error("Newton iteration failed to converge. Parent ID: $parentID.
                Time step: $n. Residual: $(norm(f(xn,system,parentID,children))).");
        end
    end

    # update junction using converged state vector
    system.branches.Q[parentID][n+2,system.solverparams.JL] = x[1];
    system.branches.A[parentID][n+2,system.solverparams.JL] = x[2];
    for i = 1:numchildren
        ID = system.branches.children[parentID][i];
        system.branches.Q[ID][n+2,1] = x[2*i+1];
        system.branches.A[ID][n+2,1] = x[2*i+2];
    end
end
