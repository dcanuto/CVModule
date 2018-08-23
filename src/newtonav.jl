function newtonav!(yout::Vector{Float64},iters::Vector{Int64},Vs::Float64,vs::Float64,
    zs::Float64,ts::Float64,tolz::Float64,state::String,Kvo::Float64,Kvc::Float64,
    Plv::Float64,Pa::Float64,Q::Float64,Vlv::Float64,zeta::Float64,W1root::Float64,W2root::Float64,
    E::Float64,V0::Float64,beta::Float64,rho::Float64,c0::Float64,A0::Float64,Aann::Float64,
    leff::Float64,Ks::Float64,th::Float64,h::Float64,f::Function,J::Function,maxiter::Int16,
    maxval::Float64,epsJ::Float64,epsN::Float64)

    # initial guess
    if state == "opening"
        zdot = (1.-zeta).*Kvo.*(Plv.-Pa.*mmHgToPa);
    else
        zdot = zeta.*Kvc.*(Plv.-Pa.*mmHgToPa);
    end
    if abs(zdot) <= tolz
        println("dζ/dt under tolerance at HR = $(60/th) bpm. tol = $tolz.
            ζ = $zeta. Switching to reduced system.")
        x0 = zeros(2);
        xx = zeros(2);
        if Q == 0
            x0[1] = -W2root + 1e-7;
        else
            x0[1] = W1root;
        end
        if Q != 0
            x0[2] = Vlv - h*Q;
        elseif state == "opening" && Q == 0
            x0[2] = Vlv;
        end

        # non-dimensionalize
        x0[1] /= vs;
        x0[2] /= Vs;
        # println(x0)
    else
        x0 = zeros(3);
        xx = zeros(3);
        if Q == 0
            x0[1] = -W2root + 1e-7;
        else
            x0[1] = W1root;
        end
        if Q != 0
            x0[2] = Vlv - h*Q;
        elseif state == "opening" && Q == 0
            x0[2] = Vlv;
        end
        if zeta == 0
            x0[3] = zs;
        elseif state == "opening"
            x0[3] = zeta + h*zdot;
        elseif state == "closing" && abs(Plv/mmHgToPa - Pa) < 20
            x0[3] = zeta + h*zdot;
        elseif state == "closing"
            x0[3] = zeta;
        end

        # non-dimensionalize
        x0[1] /= vs;
        x0[2] /= Vs;
        x0[3] /= zs;
        # println(x0)
    end

    # assign function/Jacobian handles
    # f = CVModule.fav;
    # J = CVModule.Jav;

    # setup for iterations
    xx .= x0;
    x = zeros(length(x0));
    N = 1;
    Plv = 0.;
    Pao = 0.;

    # max step size for line searches
    stpmax = 100*max(sqrt(norm(xx)),length(x0));

    # println(state)
    # Newton iterations
    while N <= maxiter
        Plv = E.*(xx[2].*Vs .- V0);
        Pao = beta.*((2.*rho./beta).*((xx[1].*vs .- W2root)./8 .+ c0).^2 .- sqrt.(A0));
        if Plv > Pao
            state = "opening";
        else
            state = "closing";
            # println(state)
        end
        # println(xx)
        # determine Jacobian, check invertibility
        JJ = J(xx,Vs,vs,ts,zs,Vlv,zeta,Q,rho,beta,W2root,c0,Kvo,Kvc,Ks,E,V0,
            A0,Aann,leff,h,state);
        D = diagm(maximum!(zeros(length(xx)),abs.(JJ)).^-1);
        # println(JJ)
        # println(D)
        # println(D*JJ)
        # println(cond(D*JJ))
        # println(det(D*JJ))
        if abs(det(D*JJ)) < epsJ
            println(JJ)
            println(D*JJ)
            print(det(D*JJ))
            error("Newton Jacobian is singular.");
        end
        # compute gradient of line search objective function
        fvec = D*f(xx,Vs,vs,ts,zs,Vlv,zeta,Q,rho,beta,W2root,c0,Kvo,Kvc,Ks,E,V0,
            A0,Aann,leff,h,state);
        # println(f(xx,system,n,state))
        # println(fvec)
        fold = 0.5*dot(fvec,fvec);
        g = (D*JJ)'*fvec;
        # compute newton step
        # println(inv(D*JJ))
        s = -inv(D*JJ)*fvec;
        # line search to update state vector
        fn, xn, check = CVModule.linesearch(xx,fold,g,s,stpmax,state,
            Vs,vs,ts,zs,Vlv,zeta,Q,rho,beta,W2root,c0,Kvo,Kvc,Ks,E,V0,
            A0,Aann,leff,h,f,J,maxiter);
        # println(xn[1]*vs)
        # check if sufficiently close to root
        JJ = J(xn,Vs,vs,ts,zs,Vlv,zeta,Q,rho,beta,W2root,c0,Kvo,Kvc,Ks,E,V0,
            A0,Aann,leff,h,state);
        D = diagm(maximum!(zeros(length(xx)),abs.(JJ)).^-1);
        fvec = D*f(xn,Vs,vs,ts,zs,Vlv,zeta,Q,rho,beta,W2root,c0,Kvo,Kvc,Ks,E,V0,
            A0,Aann,leff,h,state);
        # println(fvec)
        if norm(fvec) <= epsN*1000
            x = xn;
            x[1] = x[1]*vs;
            x[2] = x[2]*Vs;
            if length(x) == 3
                x[3] = x[3]*zs;
            end
            # println(x)
            # println(system.branches.W2root)
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
            println("Proximal Jacobian:")
            display(JJ)
            # println(D)
            # println(D*JJ)
            println("Proximal state variables (W1 (m/s), Vlv (mL), ζ):")
            println("x = $([xn[1]*vs xn[2]*Vs/cm3Tom3 xn[3]*zs])")
            println("Proximal objective function values (want zeroes): ")
            println("f = $fvec")
            println("2-norm of objective function: $(norm(fvec))")
            error("Newton iteration failed to converge.");
        end
    end

    # update based on converged solution
    yout[1] = ((2.*rho./beta).^2.*((x[1] .- W2root)./8 .+ c0).^4);
    yout[2] = 0.5.*yout[1].*(x[1] .+ W2root);
    yout[3] = x[2];
    if length(x) == 3
        if x[3] < 1e-5 && state == "closing"
            yout[4] = 0;
            println("Aortic valve closed.")
        else
            yout[4] = x[3];
        end
    else
        if zeta < 1e-6 && state == "closing"
            yout[4] = 0;
            println("Aortic valve closed.")
        else
            yout[4] = zeta;
        end
    end

    yout[5] = x[1];
end
