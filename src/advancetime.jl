function advancetime!(system::CVSystem,n::Int64;injury="no")
    while system.solverparams.numbeats < system.solverparams.numbeatstotal
        CVModule.predictorfluxes!(system,n);
        CVModule.predictorstep!(system,n);
        CVModule.correctorfluxes!(system,n);
        CVModule.correctorstep!(system,n);
        CVModule.applyendbcs!(system,n);
        CVModule.splitinvariants!(system,n);
        if injury == "no"
            CVModule.solvesplits!(system,n);
        elseif injury == "yes"
            CVModule.solvesplits!(system,n,hemoflag);
            CVModule.applytourniquet!(system,n);
        end
        CVModule.arterialpressure!(system,n);
        CVModule.regulateall!(system,n);
        # if coupleflag == "yes"
        #     CVModule.senddata(system,n,sender);
        # end
        # if n > 0 && mod(n/system.pdata.nsamp,1) == 0
        #     println(system.t[n+1])
        # end
        n+=1;
    end
    return system
end
