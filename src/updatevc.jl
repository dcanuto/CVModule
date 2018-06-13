function updatevc!(system::CVSystem,n::Int64,terms::Vector{Int64})
    ivcflow = Float64(0);
    svcflow = Float64(0);

    # sum terminal/liver outflows
    for i = 1:length(terms)
        if system.branches.group[terms[i]] == "lower"
            if terms[i] != 12 && terms[i] != 13 && terms[i] != 14 && terms[i] != 15 && terms[i] != 21 # skip portal vein, hepatic vein inputs
                ivcflow += system.branches.term[terms[i]].Q[n+1,5];
            end
        else
            svcflow += system.branches.term[terms[i]].Q[n+1,5];
        end
    end
    ivcflow += system.liver.Q[n+1,4];

    # Euler forward update to vena cavae
    system.ivc.V[n+2] = (system.ivc.V[n+1] + system.solverparams.h*(
        ivcflow - system.ivc.Q[n+1]));
    system.ivc.Q[n+2] = ((1 - system.ivc.R/system.ivc.L*
        system.solverparams.h)*system.ivc.Q[n+1]+
        system.solverparams.h/system.ivc.L*(system.ivc.P[n+1]-
        system.heart.ra.P[n+1]));
    system.ivc.P[n+2] = (1/system.ivc.C*(system.ivc.V[n+2]-
        system.ivc.V0));
    system.svc.V[n+2] = (system.svc.V[n+1] + system.solverparams.h*(
        svcflow - system.svc.Q[n+1]));
    system.svc.Q[n+2] = ((1 - system.svc.R/system.svc.L*
        system.solverparams.h)*system.svc.Q[n+1]+
        system.solverparams.h/system.svc.L*(system.svc.P[n+1]-
        system.heart.ra.P[n+1]));
    system.svc.P[n+2] = (1/system.svc.C*(system.svc.V[n+2]-
        system.svc.V0));
end
