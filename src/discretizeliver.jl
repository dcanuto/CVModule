function discretizeliver!(system::CVSystem)
    system.liver.P = zeros(system.solverparams.numsteps+1,4);
    system.liver.V = zeros(system.solverparams.numsteps+1,4);
    system.liver.Q = zeros(system.solverparams.numsteps+1,4);
end
