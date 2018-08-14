function newtondist!(yout::Vector{Float64},iters::Vector{Int64},Vs::Float64,vs::Float64,ts::Float64,
    V::Float64,C::Float64,Q::Float64,V0::Float64,beta::Float64,c0::Float64,A0::Float64,W1end::Float64,rho::Float64,
    f::Function,J::Function,maxiter::Int16,epsJ::Float64,epsN::Float64,maxval::Float64,h::Float64)

    # initial guess
    x0 = zeros(2);
    xx = zeros(2);
    x0[2] = V;
    Pest = 1/C*(x0[2] - V0);
    Aest = (Pest/beta + sqrt.(A0))^2;
    x0[1] = W1end + 8*(c0 - sqrt(0.5*beta/rho)*Aest^0.25);

    # non-dimensionalize
    x0[1] /= vs;
    x0[2] /= Vs;

    # setup for iterations
    xx .= x0;
    x = zeros(2);
    N = 1;

    # allocation
    JJ = zeros(2,2);
    D = zeros(2,2);
    d = zeros(2);
    fvec = zeros(2);
    g = zeros(2);

    # max step size for line searches
    stpmax = 100*max(sqrt(norm(xx)),2);
    # println(stpmax)

    # Newton iterations
    while N <= maxiter
        # determine Jacobian, check invertibility
        # println("Distal Jacobian:")
        J(JJ,xx,Vs,vs,ts,C,W1end,c0,rho,beta,h);
        maximum!(d,abs.(JJ))
        D[1,1] = 1/d[1];
        D[2,2] = 1/d[2];
        # println(JJ)
        # println(D)
        # println(D*JJ)
        # println(cond(D*JJ))
        # println(det(D*JJ))
        if abs(det(D*JJ)) < epsJ
            println(JJ)
            println(D*JJ)
            print(det(D*JJ))
            error("Distal Newton Jacobian is singular.");
        end
        # compute gradient of line search objective function
        # println("Objective function:")
        f(fvec,xx,D,V,Q,Vs,vs,ts,V0,rho,beta,C,W1end,c0,A0,h);
        # println(f(xx,system,n,state))
        # println(fvec)
        fold = 0.5*dot(fvec,fvec);
        g = (D*JJ)'*fvec;
        # compute newton step
        # println(inv(D*JJ))
        # println("Newton step:")
        s = -inv(D*JJ)*fvec;
        # line search to update state vector
        # println("Line search:")
        fn, xn, check = CVModule.linedist(xx,fold,g,s,stpmax,f,J,maxiter,
            V,Q,Vs,vs,ts,V0,rho,beta,C,W1end,c0,A0,h,JJ,d,fvec,D);
        # println(xn[1]*vs)
        # check if sufficiently close to root
        # println("Root checking:")
        maximum!(d,abs.(JJ))
        D[1,1] = 1/d[1];
        D[2,2] = 1/d[2];
        f(fvec,xn,D,V,Q,Vs,vs,ts,V0,rho,beta,C,W1end,c0,A0,h);
        # println(fvec)
        if norm(fvec) <= epsN
            x .= xn;
            x[1] = x[1]*vs;
            x[2] = x[2]*Vs;
            # println(inv(D*JJ))
            # println(fvec)
            iters[1] += N;
            break
        end
        if check
            test = 0.
            den = max(fn,0.5*length(xn));
            for i = 1:length(xn)
                temp = abs(g[i])*max(abs(xn[i]),1.)/den;
                if temp > test
                    test = temp;
                end
            end
            if test < 1e-6
                check = true;
                println("Warning: gradient of objective function close to zero.")
            else
                check = false;
            end
        end
        # check for divergence
        if norm(fvec,Inf) >= maxval
            iters[1] += N;
            println(xn)
            println(fvec)
            error("Newton iteration diverged.");
        end
        N+=1;
        xx = xn;
        if N == maxiter
            println(JJ)
            # println(D)
            # println(D*JJ)
            println(xn)
            println(fvec)
            println(norm(fvec))
            error("Newton iteration failed to converge.");
        end
    end

    # update based on converged solution
    # println("Update:")
    yout[1] = x[2];
    yout[2] = 1/C*(x[2] - V0);
    yout[3] = (2*rho/beta)^2*((W1end-x[1])/8 + c0)^4;
    yout[4] = yout[3]*0.5*(W1end + x[1]);
    # println("Converged distal volume, mL: $(yout[1]/cm3Tom3)")
    # println("Converged distal pressure, mmHg: $(yout[2]/mmHgToPa)")
    # println("Converged distal area, cm2: $(yout[3]/(cmTom^2))")
    # println("Converged distal flowrate, mL/s: $(yout[4]/cm3Tom3)")
end
