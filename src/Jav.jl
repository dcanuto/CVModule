function Jav(x::Vector{Float64},Vs::Float64,vs::Float64,ts::Float64,zs::Float64,
    Vlv::Float64,zeta::Float64,Q::Float64,rho::Float64,beta::Float64,W2root::Float64,
    c0::Float64,Kvo::Float64,Kvc::Float64,Ks::Float64,E::Float64,V0::Float64,
    A0::Float64,Aann::Float64,leff::Float64,h::Float64,state::String)
    # # non-dimensionalizing parameters
    # Vs = 100*cm3Tom3;
    # vs = system.branches.c0[1][end];
    # # vs = 1;
    # ts = system.heart.av.leff/vs;
    # if system.heart.av.zeta[n+1] == 0
    #     zs = 1e-5;
    # else
    #     zs = system.heart.av.zeta[n+1];
    # end

    if length(x) == 3
        J = zeros(3,3);
        J[1,1] = (-2.*ts.*vs.^5./Vs.*(rho./beta).^2.*(((x[1] .- W2root./vs)./8 .+
            c0./vs).^3.*((x[1] .+ W2root./vs)./2) .+ ((x[1] .- W2root./vs)./8 .+
            c0./vs).^4));
        J[1,2] = -1.*ts./h;
        # J11 = (-2*ts*vs^5/Vs*(system.solverparams.rho/system.branches.beta[1][end])^2*
        #     (((x[1]-system.branches.W2root/vs)/8+system.branches.c0[1][end]/vs)^3*
        #     ((x[1]+system.branches.W2root/vs)/2)+((x[1]-system.branches.W2root/vs)/8+
        #     system.branches.c0[1][end]/vs)^4));
        # J12 = -1*ts/system.solverparams.h;
        # J13 = 0;
        # J1 = [J11 J12 J13];

        if state == "opening"
            J[2,1] = (0.5.*h.*(1./zs-x[3]).*rho.*vs.^2.*Kvo.*((x[1].-W2root./vs)./8 .+ c0./vs));
            J[2,2] = Vs.*(-h.*(1./zs.-x[3]).*Kvo.*E);
            J[2,3] = (1 .+ h.*Kvo.*(E.*Vs.*(x[2] .- V0./Vs).-beta.*((2.*rho./beta).*vs.^2
                .*((x[1] .- W2root./vs)./8 .+ c0./vs).^2 .- sqrt.(A0))));
            # J21 = (0.5*system.solverparams.h*(1/zs-x[3])*system.solverparams.rho*vs^2*
            #     system.heart.av.Kvo*((x[1]-system.branches.W2root/vs)/8+
            #     system.branches.c0[1][end]/vs));
            # J22 = Vs*(-system.solverparams.h*(1/zs-x[3])*system.heart.av.Kvo*
            #     system.heart.lv.E[n+2]);
            # J23 = (1+system.solverparams.h*system.heart.av.Kvo*(system.heart.lv.E[n+2]*
            #     Vs*(x[2]-system.heart.lv.V0/Vs)-system.branches.beta[1][end]*
            #     ((2*system.solverparams.rho/system.branches.beta[1][end])*vs^2*
            #     ((x[1]-system.branches.W2root/vs)/8+system.branches.c0[1][end]/vs)^2-
            #     sqrt(system.branches.A0[1][end]))));
        else
            J[2,1] = (0.5.*h*(x[3]).*rho.*vs.^2.*Kvc.*((x[1] .- W2root./vs)./8 .+ c0./vs));
            J[2,2] = Vs.*(-h.*(x[3]).*Kvc.*E);
            J[2,3] = (1 .+ h.*Kvc.*(E.*Vs.*(x[2] .- V0./Vs) .- beta.*((2.*rho./beta).*vs.^2.*
                ((x[1].- W2root./vs)./8 .+ c0./vs).^2 .- sqrt.(A0))));
            # J21 = (0.5*system.solverparams.h*(x[3])*system.solverparams.rho*vs^2*
            #     system.heart.av.Kvc*((x[1]-system.branches.W2root/vs)/8+
            #     system.branches.c0[1][end]/vs));
            # J22 = Vs*(-system.solverparams.h*(x[3])*system.heart.av.Kvc*
            #     system.heart.lv.E[n+2]);
            # J23 = (1+system.solverparams.h*system.heart.av.Kvc*(system.heart.lv.E[n+2]*
            #     Vs*(x[2]-system.heart.lv.V0/Vs)-system.branches.beta[1][end]*
            #     ((2*system.solverparams.rho/system.branches.beta[1][end])*vs^2*
            #     ((x[1]-system.branches.W2root/vs)/8+system.branches.c0[1][end]/vs)^2-
            #     sqrt(system.branches.A0[1][end]))));
        end
        # J2 = [J21 J22 J23];

        # J31 = 0;
        if state == "closing" && x[2] > Vlv
            J[3,2]= (-x[3].*ts./h .+ Vs.*ts./(2.*Aann.*leff.*h.*zs).*
                ((Vlv./Vs .- x[2]) .- abs.(Vlv./Vs .- x[2])) .+ x[3].^2.*Ks.*E.*Vs.*ts.*zs.*
                Aann./(rho.*leff).*(Vlv./Vs .+ V0./Vs .-2.*x[2]));
            # J32 = (-x[3]*ts/system.solverparams.h+Vs*ts/(2*system.heart.av.Aann*
            #     system.heart.av.leff*system.solverparams.h*zs)*
            #     ((system.heart.lv.V[n+1]/Vs-x[2])-abs(system.heart.lv.V[n+1]/Vs-x[2]))+
            #     x[3]^2*system.heart.av.Ks*system.heart.lv.E[n+2]*Vs*ts*zs*
            #     system.heart.av.Aann/(system.solverparams.rho*system.heart.av.leff)*
            #     (system.heart.lv.V[n+1]/Vs+system.heart.lv.V0/Vs-2*x[2]));
        elseif state == "opening" || (state == "closing" && x[2] < Vlv)
            J[3,2] = (-x[3].*ts./h .+ Vs.*ts./(2.*Aann.*leff.*h.*zs).*
                ((-Vlv./Vs .+ x[2]) .- abs.(Vlv./Vs .- x[2])) .+ x[3].^2.*Ks.*E.*Vs.*ts.*zs.*
                Aann./(rho.*leff).*(Vlv./Vs .+ V0./Vs .- 2.*x[2]));
            # J32 = (-x[3]*ts/system.solverparams.h+Vs*ts/(2*system.heart.av.Aann*
            #     system.heart.av.leff*system.solverparams.h*zs)*
            #     ((-system.heart.lv.V[n+1]/Vs+x[2])-abs(system.heart.lv.V[n+1]/Vs-x[2]))+
            #     x[3]^2*system.heart.av.Ks*system.heart.lv.E[n+2]*Vs*ts*zs*
            #     system.heart.av.Aann/(system.solverparams.rho*system.heart.av.leff)*
            #     (system.heart.lv.V[n+1]/Vs+system.heart.lv.V0/Vs-2*x[2]));
        end
        if state == "opening"
            J[3,3] = ((Vlv./Vs .- x[2]).*ts./h .- Aann.*ts./(rho.*leff.*
                Kvo.*Vs).*((-zeta.*x[3].^2 .+ zs.*x[3].^3)./((1./zs .- x[3]).^2) .+
                (3.*zs.*x[3].^2 .- 2.*zeta.*x[3])./(1./zs .- x[3])) .+
                2.*Ks.*E.*Aann.*ts.*Vs.*zs./(rho.*leff).*(x[2] .- V0./Vs).*
                (Vlv./Vs .- x[2]).*x[3] .- ts./Vs.*Q);
            # J33 = ((system.heart.lv.V[n+1]/Vs-x[2])*ts/system.solverparams.h-
            #     system.heart.av.Aann*ts/(system.solverparams.rho*system.heart.av.leff*
            #     system.heart.av.Kvo*Vs)*((-system.heart.av.zeta[n+1]*x[3]^2+zs*x[3]^3)/
            #     ((1/zs-x[3])^2)+(3*zs*x[3]^2-2*system.heart.av.zeta[n+1]*x[3])/(1/zs-x[3]))+
            #     2*system.heart.av.Ks*system.heart.lv.E[n+2]*system.heart.av.Aann*ts*Vs*zs/
            #     (system.solverparams.rho*system.heart.av.leff)*
            #     (x[2]-system.heart.lv.V0/Vs)*(system.heart.lv.V[n+1]/Vs-x[2])*x[3]-
            #     ts/Vs*system.branches.Q[1][1,n+1]);
        else
            J[3,3] = ((Vlv./Vs .- x[2]).*ts./h .- Aann.*ts./(rho.*leff.*
                Kvc.*Vs).*(2.*x[3].*zs .- zeta) .+ 2.*Ks.*E.*Aann.*ts.*Vs.*zs./
                (rho.*leff).*(x[2] .- V0./Vs).*(Vlv./Vs .- x[2]).*x[3] .-
                ts./Vs.*Q);
            # J33 = ((system.heart.lv.V[n+1]/Vs-x[2])*ts/system.solverparams.h-
            #     system.heart.av.Aann*ts/(system.solverparams.rho*system.heart.av.leff*
            #     system.heart.av.Kvc*Vs)*(2*x[3]*zs-system.heart.av.zeta[n+1])+
            #     2*system.heart.av.Ks*system.heart.lv.E[n+2]*system.heart.av.Aann*ts*Vs*zs/
            #     (system.solverparams.rho*system.heart.av.leff)*
            #     (x[2]-system.heart.lv.V0/Vs)*(system.heart.lv.V[n+1]/Vs-x[2])*x[3]-
            #     ts/Vs*system.branches.Q[1][1,n+1]);
        end
        # J3 = [J31 J32 J33];

        # J = [J1;J2;J3];
    else
        J = zeros(2,2);
        J[1,1] = (-2.*ts.*vs.^5./Vs.*(rho./beta).^2.*(((x[1] .- W2root./vs)./8 .+
            c0./vs).^3.*((x[1] .+ W2root./vs)./2) .+ ((x[1] .- W2root./vs)./8 .+
            c0./vs).^4));
        # J11 = (-2*ts*vs^5/Vs*(system.solverparams.rho/system.branches.beta[1][end])^2*
        #     (((x[1]-system.branches.W2root/vs)/8+system.branches.c0[1][end]/vs)^3*
        #     ((x[1]+system.branches.W2root/vs)/2)+((x[1]-system.branches.W2root/vs)/8+
        #     system.branches.c0[1][end]/vs)^4));
        J[1,2] = -1.*ts./h;
        # J12 = -1*ts/system.solverparams.h;
        # J1 = [J11 J12];

        # J31 = 0;
        if state == "closing" && x[2] > Vlv
            J[2,2] = (-ts./h .+ Vs.*ts./(2.*Aann.*leff.*h.*zs).*((Vlv./Vs .- x[2]) .-
                abs.(Vlv./Vs .- x[2])) .+ Ks.*E.*Vs.*ts.*zs.*Aann./(rho.*leff).*
                (Vlv./Vs .+ V0./Vs .-2.*x[2]));
            # J32 = (-ts/system.solverparams.h+Vs*ts/(2*system.heart.av.Aann*
            #     system.heart.av.leff*system.solverparams.h*zs)*
            #     ((system.heart.lv.V[n+1]/Vs-x[2])-abs(system.heart.lv.V[n+1]/Vs-x[2]))+
            #     system.heart.av.Ks*system.heart.lv.E[n+2]*Vs*ts*zs*
            #     system.heart.av.Aann/(system.solverparams.rho*system.heart.av.leff)*
            #     (system.heart.lv.V[n+1]/Vs+system.heart.lv.V0/Vs-2*x[2]));
        elseif state == "opening" || (state == "closing" && x[2] < Vlv)
            J[2,2] = (-ts./h .+ Vs.*ts./(2.*Aann.*leff.*h.*zs).*((-Vlv./Vs .+ x[2]) .-
                abs.(Vlv./Vs .- x[2])) .+ Ks.*E.*Vs.*ts.*zs.*Aann./(rho.*leff).*
                (Vlv./Vs .+ V0./Vs .- 2.*x[2]));
            # J32 = (-ts/system.solverparams.h+Vs*ts/(2*system.heart.av.Aann*
            #     system.heart.av.leff*system.solverparams.h*zs)*
            #     ((-system.heart.lv.V[n+1]/Vs+x[2])-abs(system.heart.lv.V[n+1]/Vs-x[2]))+
            #     system.heart.av.Ks*system.heart.lv.E[n+2]*Vs*ts*zs*
            #     system.heart.av.Aann/(system.solverparams.rho*system.heart.av.leff)*
            #     (system.heart.lv.V[n+1]/Vs+system.heart.lv.V0/Vs-2*x[2]));
        end

        # J3 = [J31 J32];

        # J = [J1;J3];
    end
    return J
end
