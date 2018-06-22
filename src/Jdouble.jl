function Jdouble(x::Vector{Float64},beta::Vector{Float64},rho::Float64)
    J = zeros(6,6);

    J[1,1] = 1.;
    J[2,1] = rho.*x[1]./(x[2].^2);
    J[3,1] = 1.;
    J[4,1] = 1./x[2];

    J[2,2] = (0.5.*beta[1].*x[2].^-0.5 .- rho.*x[1].^2.*x[2].^-3);
    J[3,2] = J[2,2];
    J[4,2] = -x[1].*x[2].^-2 .+ sqrt.(0.5.*beta[1]./rho).*x[2].^-0.75;

    J[1,3] = -1.;
    J[2,3] = -rho.*x[3]./(x[4].^2);
    J[5,3] = 1./x[4];

    J[2,4] = (-0.5.*beta[2].*x[4].^-0.5 .+ rho.*x[3].^2.*x[4].^-3);
    J[5,4] = -x[3].*x[4].^-2 .- sqrt.(0.5.*beta[2]./rho).*x[4].^-0.75;

    J[1,5] = -1.;
    J[3,5] = -rho.*x[5]./(x[6].^2);
    J[6,5] = 1./x[6];

    J[3,6] = (-0.5.*beta[3].*x[6].^-0.5 .+ rho.*x[5].^2.*x[6].^-3);
    J[6,6] = -x[5].*x[6].^-2 .- sqrt.(0.5.*beta[3]./rho).*x[6].^-0.75;

    return J
end
