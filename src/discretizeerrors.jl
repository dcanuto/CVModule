function discretizeerrors!(system::CVSystem)
    pressure_dev = 2.5;
    area_dev = 1.36e-4;
    beta_dev = 2.11e6;
    a0_dev = 2.46e-6;
    t1_dev = 0.07;
    numarteries = 4;

    system.error.odev = [pressure_dev for i in 1:system.solverparams.JL];
    system.error.mdev = [area_dev for i in 1:system.solverparams.JL];
    for i = 1:numarteries
        push!(system.error.pdev,beta_dev);
    end
    # push!(system.error.pdev,a0_dev);
    push!(system.error.pdev,t1_dev);

    # parameter distribution smoothing
    δ = 0.985;
    system.error.a = (3*δ-1)/(2*δ);
    system.error.h = sqrt.(1-system.error.a^2);

    # parameter distribution initial mean estimates, bounds (if needed)
    for i = 1:numarteries
        push!(system.error.lb,-Inf);
        push!(system.error.ub,Inf);
        push!(system.error.pbar,2e7);
    end
    push!(system.error.lb,0.05);
    push!(system.error.ub,Inf);
    push!(system.error.pbar,0.2);

    return system
end
