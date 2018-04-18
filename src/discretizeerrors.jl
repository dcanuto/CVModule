function discretizeerrors!(system::CVSystem)
    pressure_dev = 2.5;
    area_dev = 1.36e-4;
    beta_dev = 2.11e6;
    a0_dev = 2.46e-6;
    t1_dev = 0.07;
    numarteries = 5;

    system.error.odev = [pressure_dev for i in 1:system.solverparams.JL];
    system.error.mdev = [area_dev for i in 1:system.solverparams.JL];
    # aortic β deviation
    push!(system.error.pdev,beta_dev/10);
    # brachial β deviation
    for i = 2:numarteries
        push!(system.error.pdev,beta_dev);
    end
    # push!(system.error.pdev,a0_dev);
    push!(system.error.pdev,t1_dev);

    # parameter distribution smoothing
    δ = 0.985;
    system.error.a = (3*δ-1)/(2*δ);
    system.error.h = sqrt.(1-system.error.a^2);

    # parameter distribution initial mean estimates, bounds (if needed)
    # aorta
    push!(system.error.lb,1e3);
    push!(system.error.ub,Inf);
    push!(system.error.pbar,2e6);
    # brachial artery
    for i = 2:numarteries
        push!(system.error.lb,1e3);
        push!(system.error.ub,Inf);
        push!(system.error.pbar,2e7);
    end
    # τ1
    push!(system.error.lb,0.05);
    push!(system.error.ub,Inf);
    push!(system.error.pbar,0.2);

    return system
end
