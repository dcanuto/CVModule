function advancetime!(system::CVSystem,times::CVTimer,n::Int64,nsamp::Int64,
    terms::Vector{Int64},splits::Vector{Int64};injury="no",runtype="predict")
    tic();
    if runtype == "predict"
        while system.solverparams.numbeats < system.solverparams.numbeatstotal
            if injury == "no"
                CVModule.tvdrk3!(system,times,n,splits,terms);
                tic();
                CVModule.update0d!(system,n,terms);
            elseif injury == "yes"
                CVModule.tvdrk3!(system,times,n,splits,terms,injury);
                tic();
                CVModule.update0d!(system,n,terms,injury);
            end
            times.tc += toq();
            tic();
            CVModule.arterialpressure!(system,n);
            CVModule.regulateall!(system,n,terms);
            times.tr += toq();
            if mod(n/nsamp,1) == 0
                n+=1;
                break
            end
            n+=1;
        end
    elseif runtype == "assim"
        if injury == "no"
            CVModule.tvdrk3!(system,times,n-1,splits,terms);
            tic();
            CVModule.update0d!(system,n-1,terms);
        elseif injury == "yes"
            CVModule.tvdrk3!(system,times,n-1,splits,terms,injury);
            tic();
            CVModule.update0d!(system,n-1,terms,injury);
        end
        times.tc += toq();
        tic();
        CVModule.arterialpressure!(system,n-1);
        CVModule.regulateall!(system,n-1,terms);
        times.tr += toq();
    end
    if runtype == "predict"
        times.tt += toc();
    elseif runtype == "assim"
        times.tt = toq();
    end
    return system,n,times
end
