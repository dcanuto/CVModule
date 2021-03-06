function discretizecns!(system::CVSystem)
    append!(system.cns.H,zeros(system.solverparams.numsteps+1))
    append!(system.cns.Emaxlv,zeros(system.solverparams.numsteps+1))
    append!(system.cns.Emaxrv,zeros(system.solverparams.numsteps+1))
    append!(system.cns.R2L,zeros(system.solverparams.numsteps+1))
    append!(system.cns.R3L,zeros(system.solverparams.numsteps+1))
    append!(system.cns.R2U,zeros(system.solverparams.numsteps+1))
    append!(system.cns.R3U,zeros(system.solverparams.numsteps+1))
    append!(system.cns.C4L,zeros(system.solverparams.numsteps+1))
    append!(system.cns.C5L,zeros(system.solverparams.numsteps+1))
    append!(system.cns.C4U,zeros(system.solverparams.numsteps+1))
    append!(system.cns.C5U,zeros(system.solverparams.numsteps+1))
    append!(system.cns.V4L,zeros(system.solverparams.numsteps+1))
    append!(system.cns.V5L,zeros(system.solverparams.numsteps+1))
    append!(system.cns.V4U,zeros(system.solverparams.numsteps+1))
    append!(system.cns.V5U,zeros(system.solverparams.numsteps+1))
end
