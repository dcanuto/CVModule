function Jdist!(J::Matrix{Float64},x::Vector{Float64},Vs::Float64,vs::Float64,ts::Float64,C::Float64,
    W1end::Float64,c0::Float64,rho::Float64,beta::Float64,h::Float64)
    J[1,1] = C*rho*vs^2/(2*Vs)*((W1end/vs - x[1])/8 + c0/vs)^2;
    J[1,2] = 1;
    J[2,1] = 2*ts*vs^5*rho^2/(Vs*beta^2)*(((W1end/vs - x[1])/8 + c0/vs)^3*0.5*
        (W1end/vs + x[1]) - ((W1end./vs - x[1])/8 + c0/vs)^4);
    J[2,2] = ts/h;
end
