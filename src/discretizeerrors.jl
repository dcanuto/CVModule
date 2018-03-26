function discretizeerrors!(system::CVSystem)
    pressure_dev = 2.5;
    area_dev = 1.36e-4;
    beta_dev = 1.20e8;
    a0_dev = 9.33e-5;

    system.error.odev = [pressure_dev for i in 1:system.solverparams.JL];
    system.error.mdev = [area_dev for i in 1:system.solverparams.JL];
    push!(system.error.pdev,beta_dev);
    push!(system.error.pdev,a0_dev);

    return system
end
