function linesearch(xold::Vector{Float64},fold::Float64,
    g::Vector{Float64},p::Vector{Float64},stpmax::Float64,
    state::String,Vs::Float64,vs::Float64,ts::Float64,zs::Float64,
    Vlv::Float64,zeta::Float64,Q::Float64,rho::Float64,beta::Float64,W2root::Float64,
    c0::Float64,Kvo::Float64,Kvc::Float64,Ks::Float64,E::Float64,V0::Float64,
    A0::Float64,Aann::Float64,leff::Float64,h::Float64,f::Function,J::Function,
    maxiter::Int16)
    alpha = 1e-4;
    tolx = 1e-7;
    check = false;
    # Vs = 100*cm3Tom3;
    # vs = system.branches.c0[1][end];
    # if system.heart.av.zeta[n+1] == 0
    #     zs = 1e-5;
    # else
    #     zs = system.heart.av.zeta[n+1];
    # end
    #
    # f = CVModule.fav;
    # J = CVModule.Jav;

    # scale if attempted step too large
    np = norm(p);
    # println(p)
    if np > stpmax
        p = p*stpmax/np;
    end
    # println(p)

    slope = dot(g,p);
    if slope >= 0
        print(g)
        print(p)
        error("Roundoff problem in line search.")
    end

    # compute minimum step fraction
    test = 0.;
    for i = 1:length(xold)
        temp = abs(p[i])/max(abs(xold[i]),1.);
        if temp > test
            test = temp;
        end
    end
    alamin = tolx/test;

    omega = 1;
    alam = omega; # start w/ desired fraction of Newton step
    N = 1;
    Plv = 0.;
    Pao = 0.;

    while N < maxiter
        x = xold+alam*p;
        Plv = E.*(x[2].*Vs .- V0);
        Pao = beta.*((2.*rho./beta).*((x[1].*vs .- W2root)./8 .+
            c0).^2 .- sqrt.(A0));
        if Plv > Pao
            state = "opening";
        else
            state = "closing";
            # println(state)
        end
        # println(x)
        # println(xold)
        # println(alam)
        # println(p)
        JJ = J(x,Vs,vs,ts,zs,Vlv,zeta,Q,rho,beta,W2root,c0,Kvo,Kvc,Ks,E,V0,
            A0,Aann,leff,h,state);
        # println(JJ)
        D = diagm(maximum!(zeros(length(x)),abs.(JJ)).^-1);
        fvec = D*f(x,Vs,vs,ts,zs,Vlv,zeta,Q,rho,beta,W2root,c0,Kvo,Kvc,Ks,E,V0,
            A0,Aann,leff,h,state);
        # println(D)
        # println(fvec)
        fn = 0.5*dot(fvec,fvec);
        if alam < alamin
            x = xold;
            # println(alam)
            # println(alamin)
            # println(fn)
            # println(fold+alpha*alam*slope)
            println("Δx converged in proximal line search. Verify in Newton loop.")
            check = true;
            return fn,x,check
        elseif fn < fold+alpha*alam*slope
            # println(fn)
            # println(fold)
            # println(alpha*alam*slope)
            return fn,x,check
        else
            if alam == omega
                tmplam = -slope/(2*(fn-fold-slope));
            else
                rhs1 = fn-fold-alam*slope
                rhs2 = fn2-fold-alam2*slope
                a = (rhs1/alam^2-rhs2/alam2^2)/(alam-alam2);
                b = (-alam2*rhs1/alam^2+alam*rhs2/alam2^2)/(alam-alam2);
                if a == 0
                    tmplam = -slope/(2*b);
                else
                    disc = b^2-3*a*slope;
                    if disc < 0
                        tmplam = 0.5*alam;
                    elseif b <= 0
                        tmplam = (-b+sqrt(disc)/(3*a));
                    else
                        tmplam = -slope/(b+sqrt(disc));
                    end
                end
                if tmplam > 0.5*alam
                    tmplam = 0.5*alam;
                end
            end
        end
        alam2 = alam;
        fn2 = fn;
        alam = max(tmplam,0.1*alam);
        N+=1;
        if N == maxiter
            println(xn)
            println(f(xn,Vs,vs,ts,zs,Vlv,zeta,Q,rho,beta,W2root,c0,Kvo,Kvc,Ks,E,V0,
                A0,Aann,leff,h,state))
            error("Proximal line search iteration failed to converge.");
        end
    end
end
