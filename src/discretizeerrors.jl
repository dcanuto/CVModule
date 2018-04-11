function discretizeerrors!(system::CVSystem)
    pressure_dev = 2.5;
    area_dev = 1.36e-4;
    beta_dev = 2.11e6;
    a0_dev = 2.46e-6;
    Elv_dev = 4e7;
    Erv_dev = 5.33e6;
    t1_dev = 0.07

    system.error.odev = [pressure_dev for i in 1:system.solverparams.JL];
    system.error.mdev = [area_dev for i in 1:system.solverparams.JL];
    push!(system.error.pdev,beta_dev);
    # push!(system.error.pdev,a0_dev);
    # push!(system.error.pdev,Elv_dev);
    # push!(system.error.pdev,Erv_dev);
    push!(system.error.pdev,t1_dev);

    return system
end
