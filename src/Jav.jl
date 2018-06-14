function Jav(x::Vector{Float64},Vs::Float64,vs::Float64,ts::Float64,zs::Float64,
    Vlv::Float64,zeta::Float64,Q::Float64,rho::Float64,beta::Float64,W2root::Float64,
    c0::Float64,Kvo::Float64,Kvc::Float64,Ks::Float64,E::Float64,V0::Float64,
    A0::Float64,Aann::Float64,leff::Float64,h::Float64,state::String)

    if length(x) == 3
        J = zeros(3,3);
        J[1,1] = (-2.*ts.*vs.^5./Vs.*(rho./beta).^2.*(((x[1] .- W2root./vs)./8 .+
            c0./vs).^3.*((x[1] .+ W2root./vs)./2) .+ ((x[1] .- W2root./vs)./8 .+
            c0./vs).^4));
        J[1,2] = -1.*ts./h;

        if state == "opening"
            J[2,1] = (0.5.*h.*(1./zs-x[3]).*rho.*vs.^2.*Kvo.*((x[1].-W2root./vs)./8 .+ c0./vs));
            J[2,2] = Vs.*(-h.*(1./zs.-x[3]).*Kvo.*E);
            J[2,3] = (1 .+ h.*Kvo.*(E.*Vs.*(x[2] .- V0./Vs).-beta.*((2.*rho./beta).*vs.^2
                .*((x[1] .- W2root./vs)./8 .+ c0./vs).^2 .- sqrt.(A0))));
        else
            J[2,1] = (0.5.*h*(x[3]).*rho.*vs.^2.*Kvc.*((x[1] .- W2root./vs)./8 .+ c0./vs));
            J[2,2] = Vs.*(-h.*(x[3]).*Kvc.*E);
            J[2,3] = (1 .+ h.*Kvc.*(E.*Vs.*(x[2] .- V0./Vs) .- beta.*((2.*rho./beta).*vs.^2.*
                ((x[1].- W2root./vs)./8 .+ c0./vs).^2 .- sqrt.(A0))));
        end
        if state == "closing" && x[2] > Vlv # closing with backflow
            J[3,2]= (-x[3].*ts./h .+ Vs.*ts./(2.*Aann.*leff.*h.*zs).*
                ((Vlv./Vs .- x[2]) .- abs.(Vlv./Vs .- x[2])) .+ x[3].^2.*Ks.*E.*Vs.*ts.*zs.*
                Aann./(rho.*leff).*(Vlv./Vs .+ V0./Vs .-2.*x[2]));
        elseif state == "opening" || (state == "closing" && x[2] < Vlv) # closing w/ forward flow
            J[3,2] = (-x[3].*ts./h .+ Vs.*ts./(2.*Aann.*leff.*h.*zs).*
                ((-Vlv./Vs .+ x[2]) .- abs.(Vlv./Vs .- x[2])) .+ x[3].^2.*Ks.*E.*Vs.*ts.*zs.*
                Aann./(rho.*leff).*(Vlv./Vs .+ V0./Vs .- 2.*x[2]));
        end
        if state == "opening"
            J[3,3] = ((Vlv./Vs .- x[2]).*ts./h .- Aann.*ts./(rho.*leff.*
                Kvo.*Vs).*((-zeta.*x[3].^2 .+ zs.*x[3].^3)./((1./zs .- x[3]).^2) .+
                (3.*zs.*x[3].^2 .- 2.*zeta.*x[3])./(1./zs .- x[3])) .+
                2.*Ks.*E.*Aann.*ts.*Vs.*zs./(rho.*leff).*(x[2] .- V0./Vs).*
                (Vlv./Vs .- x[2]).*x[3] .- ts./Vs.*Q);
        else
            J[3,3] = ((Vlv./Vs .- x[2]).*ts./h .- Aann.*ts./(rho.*leff.*
                Kvc.*Vs).*(2.*x[3].*zs .- zeta) .+ 2.*Ks.*E.*Aann.*ts.*Vs.*zs./
                (rho.*leff).*(x[2] .- V0./Vs).*(Vlv./Vs .- x[2]).*x[3] .-
                ts./Vs.*Q);
        end
    else # reduced system (constant aortic valve state)
        J = zeros(2,2);
        J[1,1] = (-2.*ts.*vs.^5./Vs.*(rho./beta).^2.*(((x[1] .- W2root./vs)./8 .+
            c0./vs).^3.*((x[1] .+ W2root./vs)./2) .+ ((x[1] .- W2root./vs)./8 .+
            c0./vs).^4));
        J[1,2] = -1.*ts./h;

        if state == "closing" && x[2] > Vlv # closing w/ backflow
            J[2,2] = (-ts./h .+ Vs.*ts./(2.*Aann.*leff.*h.*zs).*((Vlv./Vs .- x[2]) .-
                abs.(Vlv./Vs .- x[2])) .+ Ks.*E.*Vs.*ts.*zs.*Aann./(rho.*leff).*
                (Vlv./Vs .+ V0./Vs .-2.*x[2])); # closing w/ forward flow
        elseif state == "opening" || (state == "closing" && x[2] < Vlv)
            J[2,2] = (-ts./h .+ Vs.*ts./(2.*Aann.*leff.*h.*zs).*((-Vlv./Vs .+ x[2]) .-
                abs.(Vlv./Vs .- x[2])) .+ Ks.*E.*Vs.*ts.*zs.*Aann./(rho.*leff).*
                (Vlv./Vs .+ V0./Vs .- 2.*x[2]));
        end
    end
    return J
end
