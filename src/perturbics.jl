function perturbics!(system::CVSystem)
    # nptbs = 100; # number of perturbations
    # fracmax = 0.1; # max. fraction of compartmental blood volume to switch
    # termidxs = []; # indices of branches without children
    #
    # for i = 1:length(system.branches.ID)
    #     if isempty(system.branches.children[i])
    #         push!(termidxs,i);
    #     end
    # end
    #
    # for i = 1:nptbs
    #     ptbidxs = shuffle(termidxs)[1:2];
    #     frac = rand()*fracmax;
    #     for j = 2:length(system.branches.term[ptbidxs[1]].V[1,:])
    #         system.branches.term[ptbidxs[2]].V[1,j] +=
    #             frac*system.branches.term[ptbidxs[1]].V[1,j];
    #         system.branches.term[ptbidxs[1]].V[1,j] *= (1-frac);
    #         system.branches.term[ptbidxs[1]].P[1,j] =
    #             1/system.branches.term[ptbidxs[1]].C[j]*
    #             (system.branches.term[ptbidxs[1]].V[1,j]-
    #             system.branches.term[ptbidxs[1]].V0[1]);
    #         system.branches.term[ptbidxs[2]].P[1,j] =
    #             1/system.branches.term[ptbidxs[2]].C[j]*
    #             (system.branches.term[ptbidxs[2]].V[1,j]-
    #             system.branches.term[ptbidxs[2]].V0[1]);
    #     end
    # end
    # system.branches.beta[61][end] = rand(Normal(2.11e7,system.error.pdev[1]));
    # system.branches.A0[61][end] = rand(Normal(2.46e-5,system.error.pdev[2]));

    # system.branches.beta[61][end] = rand(Normal(4e7,system.error.pdev[1]));
    system.branches.A0[61][end] = rand(Normal(1e-5,system.error.pdev[2]));

    # system.branches.beta[61][end] = rand(Truncated(Normal(9.13e7,system.error.pdev[1]),2.33e6,4.89e8));
    # system.branches.A0[61][end] = rand(Truncated(Normal(4.44e-5,system.error.pdev[2]),2.01e-6,6.61e-4));
    system.branches.A[61][1,:] = (70*mmHgToPa/system.branches.beta[61][end] +
        sqrt(system.branches.A0[61][end]))^2;

    return system
end
