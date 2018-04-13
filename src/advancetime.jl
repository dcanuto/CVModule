function advancetime!(system::CVSystem,n::Int64;injury="no",runtype="predict")
    if runtype == "predict" # normal time integration
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
            if mod(n/system.pdata.nsamp,1) == 0
                n+=1
                break
            end
            n+=1;
        end
    elseif runtype == "assim" # forecast from analyzed parameters
        CVModule.predictorfluxes!(system,n-1);
        CVModule.predictorstep!(system,n-1);
        CVModule.correctorfluxes!(system,n-1);
        CVModule.correctorstep!(system,n-1);
        CVModule.applyendbcs!(system,n-1);
        CVModule.splitinvariants!(system,n-1);
        if injury == "no"
            CVModule.solvesplits!(system,n-1);
        elseif injury == "yes"
            CVModule.solvesplits!(system,n-1,hemoflag);
            CVModule.applytourniquet!(system,n-1);
        end
        CVModule.arterialpressure!(system,n-1);
        CVModule.regulateall!(system,n-1);
    end
    return system, n
end
