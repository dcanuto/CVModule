function newton!(iters::Vector{Int64},Q::Vector{Float64},A::Vector{Float64},
    children::Vector{Int64},Qold::Vector{Float64},Aold::Vector{Float64},W::Vector{Float64},
    beta::Vector{Float64},A0::Vector{Float64},c0::Vector{Float64},rho::Float64,
    f::Function,J::Function,maxiter::Int16,maxval::Float64,epsJ::Float64,epsN::Float64)

    numchildren = length(children);

    # initial guess state vector (flow rate/area pairs)
    x0 = zeros(2*numchildren+2);
    xx = zeros(2*numchildren+2);
    x0[1] = Qold[1];
    x0[2] = Aold[1];
    for i = 1:numchildren
        x0[2*i+1] = Qold[i+1];
        x0[2*i+2] = Aold[i+1];
    end

    xx .= x0;
    x = zeros(length(x0));
    N = 1; # iteration counter

    # Newton iterations
    while N <= maxiter
        # determine Jacobian, check invertibility
        JJ = J(xx,beta,rho)
        if abs(det(JJ)) < epsJ
            error("Newton Jacobian is singular.");
        end
        # update state vector
        xn = xx - inv(JJ)*f(xx,beta,A0,rho,c0,W);
        # check if sufficiently close to root
        if norm(f(xn,beta,A0,rho,c0,W)) <= epsN
            x = xn;
            # println(x)
            iters[1] += N;
            break
        end
        # check for divergence
        if norm(f(xn,beta,A0,rho,c0,W)) >= maxval
            iters[1] += N;
            error("Newton iteration diverged.");
        end
        N+=1;
        xx = xn;
        if N == maxiter
            error("Newton iteration failed to converge.");
        end
    end

    # update junction using converged state vector
    Q[1] = x[1];
    A[1] = x[2];
    for i = 1:numchildren
        Q[i+1] = x[2*i+1];
        A[i+1] = x[2*i+2];
    end
end
