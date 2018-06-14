function fav(x::Vector{Float64},Vs::Float64,vs::Float64,ts::Float64,zs::Float64,
    Vlv::Float64,zeta::Float64,Q::Float64,rho::Float64,beta::Float64,W2root::Float64,
    c0::Float64,Kvo::Float64,Kvc::Float64,Ks::Float64,E::Float64,V0::Float64,
    A0::Float64,Aann::Float64,leff::Float64,h::Float64,state::String)

    if length(x) == 3
        f = zeros(3);
        f[1] = ((Vlv./Vs .- x[2]).*ts./h .- ts.*vs.^5./Vs.*(2.*rho./beta).^2.*
            ((x[1] .- W2root./vs)./8 .+ c0./vs).^4.*((x[1] .+ W2root./vs)./2));
        if state == "opening"
            f[2] = (x[3] .- h.*((1./zs .- x[3]).*Kvo.*(E.*Vs.*(x[2] .- V0./Vs) .-
                beta.*((2.*rho./beta).*vs.^2.*((x[1] .- W2root./vs)./8 .+
                c0./vs).^2 .- sqrt.(A0)))) .- zeta./zs);
            f[3] = (x[3].*(Vlv./Vs .- x[2]).*ts./h .- Aann.*ts./(rho.*leff.*Vs).*
                ((zs.*x[3].^3 .- zeta.*(x[3]).^2)./(Kvo.*(1./zs .- x[3])) .- 1./zs.*
                rho.*Vs.^2./(2.*Aann.^2.*h).*(Vlv./Vs .- x[2]).*
                abs.(Vlv./Vs .- x[2]) .- zs.*x[3].^2.*Ks.*E.*Vs.^2.*(x[2] .- V0./Vs).*
                (Vlv./Vs .- x[2])) .- ts./Vs.*x[3].*Q);
        else
            f[2] = (x[3] .- h.*((x[3]).*Kvc.*(E.*Vs.*(x[2] .- V0./Vs) .-
                beta.*((2.*rho./beta).*vs.^2.*((x[1] .- W2root./vs)./8 .+
                c0./vs).^2 .- sqrt.(A0)))) .- zeta./zs);
            f[3] = (x[3].*(Vlv./Vs .- x[2]).*ts./h .- Aann.*ts./(rho.*leff.*Vs).*
                ((zs.*x[3].^2 .- zeta.*x[3])./Kvc .- 1./zs.*rho.*Vs.^2./(2.*Aann.^2.*
                h).*(Vlv./Vs .- x[2]).*abs.(Vlv./Vs .- x[2]) .- zs.*x[3].^2.*Ks.*
                E.*Vs.^2.*(x[2] .- V0./Vs).*(Vlv./Vs .- x[2])) .- ts./Vs.*x[3].*Q);
        end
    else
        f = zeros(2);
        f[1] = ((Vlv./Vs .- x[2]).*ts./h .- ts.*vs.^5./Vs.*(2.*rho./beta).^2.*
            ((x[1] .- W2root./vs)./8 .+ c0./vs).^4.*((x[1] .+ W2root./vs)./2));
        f[2] = ((Vlv./Vs .- x[2]).*ts./h .- Aann.*ts./(rho.*leff.*Vs).*
            (-1./zs.*rho.*Vs.^2./(2.*Aann.^2.*h).*(Vlv./Vs .- x[2]).*
            abs.(Vlv./Vs .- x[2]) .- zs.*Ks.*E.*Vs.^2.*(x[2] .- V0./Vs).*
            (Vlv./Vs .- x[2])) .- ts./Vs.*Q);
    end
    return f
end
