function splitinvariants!(W::Vector{Float64},Q::Vector{Float64},A::Vector{Float64},
    Qo::Vector{Float64},Ao::Vector{Float64},beta::Vector{Float64},c0::Vector{Float64},
    rho::Float64,k::Vector{Float64},h::Float64,diff::Float64,mu::Float64,children::Vector{Int64})
    lf = Float64;
    lb = Float64;
    W1 = Vector{Float64}[];
    dW1 = Float64;
    W2 = Vector{Float64}[];
    dW2 = Float64;

    # forward wave speed for parent's distal end
    lf = (Q[1]./A[1] .+ sqrt.(0.5.*beta[1]./rho).*A[1].^0.25);
    # parent's right running invariant at current time
    W1 = ([Qo[1]./Ao[1] .+ 4.*(sqrt.(0.5.*beta[1]./rho).*Ao[1]^0.25 .- c0[1])]);
    push!(W1,Q[1]./A[1] .+ 4.*(sqrt.(0.5.*beta[1]./rho).*A[1].^0.25 .- c0[1]));
    dW1 = (W1[2] - W1[1])/k[1];
    dQ = (Q[1] .- Qo[1])./k[1];
    dA = (A[1] .- Ao[1])./k[1];
    # advance parent's distal boundary W1 by linear interp.
    Qint = Q[1] .- dQ.*lf.*h;
    Aint = A[1] .- dA.*lf.*h;
    W[1] = W1[2] .- dW1.*lf.*h .- diff.*pi.*mu./rho.*Qint./Aint.^2.*h;
    for j = 2:(length(children)+1)
        lb = (Q[j]./A[j] .- sqrt.(0.5.*beta[j]./rho).*A[j].^0.25);
        W2 = ([Q[j]./A[j] .- 4.*(sqrt.(0.5.*beta[j]./rho).*A[j].^0.25 .- c0[j])]);
        push!(W2,Qo[j]./Ao[j] .- 4.*(sqrt.(0.5.*beta[j]./rho).*Ao[j].^0.25 .- c0[j]))
        dW2 = (W2[2] .- W2[1])./k[j];
        dQ = (Qo[j] .- Q[j])./k[j];
        dA = (Ao[j] .- A[j])./k[j];
        # advance child's proximal boundary invariant
        Qint = Q[j] .- dQ.*lb.*h;
        Aint = A[j] .- dA.*lb.*h;
        W[j] = W2[1] .- dW2.*lb.*h .- diff.*pi.*mu./rho.*Qint./Aint.^2.*h;
    end
end
