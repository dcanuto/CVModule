function fdist(x::Vector{Float64},V::Float64,Q::Float64,Vs::Float64,vs::Float64,
    ts::Float64,V0::Float64,rho::Float64,beta::Float64,C::Float64,W1end::Float64,
    c0::Float64,A0::Float64,h::Float64)

    f = zeros(2);

    f[1] = x[2] .- V0./Vs .- beta.*C./Vs.*(vs.^2.*(2.*rho./beta).*((W1end./vs .- x[1])./8 .+
        c0./vs).^2 .- sqrt.(A0));
    f[2] = (x[2] .- V./Vs)./h.*ts .- ts./Vs.*(vs.^5.*(2.*rho./beta).^2.*((W1end./vs .- x[1])./8 .+
        c0./vs).^4.*0.5.*(W1end./vs .+ x[1])) .+ ts/Vs.*Q;

    return f
end
