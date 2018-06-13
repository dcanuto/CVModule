function rootinvariant!(W2root::Vector{Float64},Q::Float64,A::Float64,Qp1::Float64,Ap1::Float64,
    c0::Float64,beta::Float64,rho::Float64,mu::Float64,diff::Float64,k::Float64,h::Float64)
    lb = Float64;
    W2 = Vector{Float64}[];
    dW2 = Float64;
    dQ = Float64;
    dA = Float64;
    Qint = Float64;
    Aint = Float64;

    # backward wave speed at proximal end
    lb = 0.5.*((Qp1./Ap1 .- sqrt.(0.5.*beta./rho).*Ap1.^0.25) .+
        (Q./A .- sqrt.(0.5.*beta./rho).*A.^0.25));
    # left-running invariant at current time
    W2 = [Q./A .- 4.*(sqrt.(0.5.*beta./rho).*A.^0.25 .- c0)];
    push!(W2,Qp1./Ap1 .- 4.*(sqrt.(0.5.*beta./rho).*Ap1.^0.25 .- c0))
    dW2 = (W2[2] .- W2[1])./k;
    # change in solution variables
    dQ = (Qp1 .- Q)./k;
    dA = (Ap1 .- A)./k;
    # interpolated state variables
    Qint = Q .- dQ.*lb.*h;
    Aint = A .- dA.*lb.*h;
    # advance proximal boundary invariant
    W2root[1] = W2[1] .- dW2.*lb.*h .- diff.*mu./rho.*h.*(Qint./Aint.^2);
end
