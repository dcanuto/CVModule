function matchpdata!(system::CVSystem)
    print("Old Δt: ")
    println(system.solverparams.h)
    system.pdata.nsamp = Int(cld(system.pdata.dtsamp,system.solverparams.h)) + 1;
    print("Number of time steps between samples: ")
    println(system.pdata.nsamp)
    system.solverparams.h = system.pdata.dtsamp/system.pdata.nsamp;
    print("New Δt matched to data sampling rate: ")
    println(system.solverparams.h)
    print("Sampling interval: ")
    println(system.pdata.dtsamp)
    print("Steps between samples × New Δt (should equal sampling interval): ")
    println(system.pdata.nsamp*system.solverparams.h)
end
