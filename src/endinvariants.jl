function endinvariants!(W1end::Vector{Float64},Q::Float64,A::Float64,Qm1::Float64,Am1::Float64,
    c0::Float64,beta::Float64,rho::Float64,mu::Float64,diff::Float64,k::Float64,h::Float64)
    lf = Float64;
    W1 = Vector{Float64}[];
    dW1 = Float64;
    dQ = Float64;
    dA = Float64;
    Qint = Float64;
    Aint = Float64;

    # forward wave speed
    lf = (Q./A .+ sqrt.(0.5.*beta./rho).*A.^0.25);
    # right-running invariants at current time
    W1 = [Qm1./Am1 .+ 4.*(sqrt.(0.5.*beta./rho).*Am1.^0.25 .- c0)];
    push!(W1,Q./A .+ 4.*(sqrt.(0.5.*beta./rho).*A.^0.25 .- c0));
    dW1 = (W1[2] .- W1[1])./k;
    # change in solution variables
    dQ = (Qm1 .- Q)./k;
    dA = (Am1 .- A)./k;
    # interpolated solution variables
    Qint = Qm1 .- dQ.*lf.*h;
    Aint = Am1 .- dA.*lf.*h;
    # update right-running invariant by linear interp.
    W1end[1] = W1[2] .- dW1.*lf.*h .- diff.*mu./rho.*h.*(Qint./Aint.^2);
end
