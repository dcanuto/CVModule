function coupleproximal!(system::CVSystem,n::Int64,h::Float64)
    if (system.heart.lv.P[n+1]/mmHgToPa < system.branches.P[1][1,n+1] &&
        system.heart.av.zeta[n+1] == 0)
        AV = "closed";
    elseif (system.heart.av.zeta[n+1] >= 0 && system.heart.lv.P[n+1]/mmHgToPa >
        system.branches.P[1][1,n+1])
        AV = "opening";
    elseif system.heart.lv.P[n+1]/mmHgToPa < system.branches.P[1][1,n+1]
        AV = "closing";
    end

    if AV == "closed"
        # # update left-running invariant by extrapolation
        # ret1 = Float64[0];
        # CVModule.rootinvariant!(ret1,system.branches.Q[1][1,n+1],system.branches.A[1][1,n+1],
        #     system.branches.Q[1][2,n+1],system.branches.A[1][2,n+1],system.branches.c0[1][end],
        #     system.branches.beta[1][end],system.solverparams.rho,system.solverparams.mu,
        #     system.solverparams.diffusioncoeff,system.branches.k[1],system.solverparams.h);
        # system.branches.W2root = ret1[1];
        # update right-running invariant
        system.branches.W1root = -system.branches.W2root;
        # update proximal A, Q w/ invariants
        system.branches.A[1][1,n+2] = ((2*system.solverparams.rho/
            system.branches.beta[1][end])^2*(0.125*(system.branches.W1root-
            system.branches.W2root) + system.branches.c0[1][end])^4);
        system.branches.Q[1][1,n+2] = (system.branches.A[1][1,n+2]*0.5*(
            system.branches.W1root + system.branches.W2root));
        system.heart.lv.V[n+2] = (system.heart.lv.V[n+1]+
            system.solverparams.h*(system.heart.la.Q[n+1]-
            system.branches.Q[1][1,n+2]));
    else
        # ventricular elastance at next time step
        CVModule.elastancefn!(system,n+1,h);
        # # update W2 by extrapolation
        # ret1 = Float64[0];
        # CVModule.rootinvariant!(ret1,system.branches.Q[1][1,n+1],system.branches.A[1][1,n+1],
        #     system.branches.Q[1][2,n+1],system.branches.A[1][2,n+1],system.branches.c0[1][end],
        #     system.branches.beta[1][end],system.solverparams.rho,system.solverparams.mu,
        #     system.solverparams.diffusioncoeff,system.branches.k[1],system.solverparams.h);
        # system.branches.W2root = ret1[1];
        # non-dimensionalizing parameters
        Vs = 100*cm3Tom3;
        vs = system.branches.c0[1][end];
        # ts = system.heart.av.leff/vs;
        ts = system.branches.k[1]/vs;

        if system.heart.av.zeta[n+1] == 0
            zs = 1e-5;
        else
            zs = system.heart.av.zeta[n+1];
        end
        # tolerance on dζ/dt (needs to be lower for higher HR)
        tolz = 1e-12 + (0.8 - system.heart.activation.th[end])*2.5e-13;
        # Newton iterations
        yout = zeros(5);
        iters = Int64[0];
        CVModule.newtonav!(yout,iters,Vs,vs,zs,ts,tolz,AV,system.heart.av.Kvo,system.heart.av.Kvc,
            system.heart.lv.P[n+1],system.branches.P[1][1,n+1],system.branches.Q[1][1,n+1],
            system.heart.lv.V[n+1],system.heart.av.zeta[n+1],system.branches.W1root,
            system.branches.W2root,system.heart.lv.E[n+2],system.heart.lv.V0,
            system.branches.beta[1][end],system.solverparams.rho,system.branches.c0[1][end],
            system.branches.A0[1][end],system.heart.av.Aann,system.heart.av.leff,
            system.heart.av.Ks,system.heart.activation.th[end],system.solverparams.h,CVModule.fav,CVModule.Jav,
            system.solverparams.maxiter,system.solverparams.maxval,
            system.solverparams.epsJ,system.solverparams.epsN)
        system.solverparams.totaliter += iters[1];
        system.branches.A[1][1,n+2] = yout[1];
        system.branches.Q[1][1,n+2] = yout[2];
        system.heart.lv.V[n+2] = yout[3];
        system.heart.av.zeta[n+2] = yout[4];
        system.branches.W1root = yout[5];
    end
end
