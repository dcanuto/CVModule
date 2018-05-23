function updateliver!(system::CVSystem,n::Int64)
    inflow = Float64(0);

    # sum liver inflows
    for i = 1:length(system.branches.ID)
        if isempty(system.branches.children[i])
            if i == 12 || i == 13 || i == 15 # portal vein inputs
                inflow += system.branches.term[i].Q[n+1,5];
            end
        end
    end

    # Euler forward update to liver
    system.liver.V[n+2,1] = (system.liver.V[n+1,1] + system.solverparams.h*(
        inflow - system.liver.Q[n+1,1]));
    system.liver.Q[n+2,1] = ((1 - system.liver.R[1]/system.liver.L*
        system.solverparams.h)*system.liver.Q[n+1,1]+
        system.solverparams.h/system.liver.L*(system.liver.P[n+1,1]-
        system.liver.P[n+1,2]));
    system.liver.V[n+2,2] = (system.liver.V[n+1,2]+
        system.solverparams.h*(system.liver.Q[n+1,1] - system.liver.Q[n+1,2]));
    system.liver.V[n+2,3] = (system.liver.V[n+1,3]+
        system.solverparams.h*(system.liver.Q[n+1,2] - system.liver.Q[n+1,3]));
    system.liver.V[n+2,4] = (system.liver.V[n+1,4]+
        system.solverparams.h*(system.liver.Q[n+1,3] - system.liver.Q[n+1,4]));

    # note pressure, flowrate updated algebraically except portal vein compartment
    system.liver.P[n+2,1] = (1/system.liver.C[1]*
        (system.liver.V[n+2,1] - system.liver.V0[1]));
    system.liver.P[n+2,2] = (1/system.liver.C[2]*
        (system.liver.V[n+2,2] - system.liver.V0[2]));
    system.liver.P[n+2,3] = (1/system.liver.C[3]*
        (system.liver.V[n+2,3] - system.liver.V0[3]));
    system.liver.P[n+2,4] = (1/system.liver.C[4]*
        (system.liver.V[n+2,4] - system.liver.V0[4]));
    system.liver.Q[n+2,2] = ((system.liver.P[n+2,2] - system.liver.P[n+2,3])/
        system.liver.R[2]);
    system.liver.Q[n+2,3] = ((system.liver.P[n+2,3] - system.liver.P[n+2,4])/
        system.liver.R[3]);
    system.liver.Q[n+2,4] = ((system.liver.P[n+2,4] - system.ivc.P[n+2])/
        system.liver.R[4]);
end
