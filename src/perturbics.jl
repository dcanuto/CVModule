function perturbics!(system::CVSystem)
    for i = 61:64
        # perturb brahcial artery beta, A0
        system.branches.beta[i][end] = CVModule.paramwalk(system,system.branches.beta[i][end],system.error.pbar[i-60]);
        # system.branches.beta[i][end] = rand(Normal(2e7,system.error.pdev[i-60]));
        # system.branches.A0[i][end] = rand(Normal(1e-5,system.error.pdev[2]));

        # set area to give 70 mmHg initial pressure
        # system.branches.A[i][1,:] = (70*mmHgToPa/system.branches.beta[i][end] +
        #     sqrt(system.branches.A0[i][end]))^2;
        system.branches.A[i][1,:] = (50*mmHgToPa/system.branches.beta[i][end] +
            sqrt(system.branches.A0[i][end]))^2;
    end

    # # perturb ventricular elastance
    # system.heart.lv.E[1] = rand(Normal(system.heart.lv.E[1],0.1*system.heart.lv.E[1]));
    # system.heart.rv.E[1] = rand(Normal(system.heart.rv.E[1],0.1*system.heart.rv.E[1]));
    #
    # # consistent ventricular volume (assuming initial P is roughly correct)
    # system.heart.lv.V[1] = system.heart.lv.P[1]/system.heart.lv.E[1] + system.heart.lv.V0;
    # system.heart.rv.V[1] = system.heart.rv.P[1]/system.heart.rv.E[1] + system.heart.rv.V0;
    system.heart.activation.tau1 = rand(Normal(0.2,system.error.pdev[end]));
    t = linspace(0,system.heart.activation.th[end],10000);
    g1 = (t/system.heart.activation.tau1).^system.heart.activation.m1;
    g2 = (t/system.heart.activation.tau2).^system.heart.activation.m2;
    system.heart.activation.k[end] = maximum((g1./(1+g1)).*(1./(1+g2)))^-1;

    return system
end
