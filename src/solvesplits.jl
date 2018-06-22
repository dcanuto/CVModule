function solvesplits!(iters::Vector{Int64},Q::Vector{Float64},A::Vector{Float64},
    children::Vector{Int64},Qold::Vector{Float64},Aold::Vector{Float64},W::Vector{Float64},
    beta::Vector{Float64},A0::Vector{Float64},c0::Vector{Float64},rho::Float64,
    f::Function,J::Function,maxiter::Int16,maxval::Float64,epsJ::Float64,epsN::Float64,hemoflag="no",
    injured = false,h = 0.,Ph = 0.,Vloss = [0.])
    if hemoflag == "no"
        CVModule.newton!(iters,Q,A,children,Qold,Aold,W,beta,A0,c0,rho,f,J,
            maxiter,maxval,epsJ,epsN);
    elseif hemoflag == "yes"
        if injured == false
            CVModule.newton!(iters,Q,A,children,Qold,Aold,W,beta,A0,c0,rho,f,J,
                maxiter,maxval,epsJ,epsN);
        elseif injured == true
            yout = zeros(3*(length(children)+1)+1);
            CVModule.model1dhemo!(yout,children,Ph,beta,rho,c0,A0,h,W[1],W[2:end])
            A[1] = yout[1];
            W[1] = yout[2];
            Q[1] = yout[3];
            Vloss[1] = yout[4];
            for i = 1:length(children)
                A[i+1] = yout[5+3*(i-1)];
                W[i+1] = yout[6+3*(i-1)];
                Q[i+1] = yout[7+3*(i-1)];
            end
        end
    end
end
