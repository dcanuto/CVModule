importall CVModule

type CVTimer
    tfd::Float64
    tc::Float64
    ts::Float64
    tr::Float64
    tt::Float64
    function CVTimer()
        this = new();
        tfd = 0;
        tc = 0;
        ts = 0;
        tr = 0;
        tt = 0;
        return this
    end
end

function main()

# filename = "arterytree.csv";j
filename = "test2.mat";
rstflag = "yes"
hemoflag = "no"
saveflag = "yes"
coupleflag = "no"

system = CVModule.buildall(filename;numbeatstotal=1,restart=rstflag,injury=hemoflag);

times = CVTimer();

n = system.solverparams.nstart;

if coupleflag == "yes"
    ctx=ZMQ.Context();
    sender=ZMQ.Socket(ctx,ZMQ.REQ);
    ZMQ.connect(sender,"tcp://127.0.0.1:5555");
end

tic();
while system.solverparams.numbeats < system.solverparams.numbeatstotal
    tic();
    CVModule.predictorfluxes!(system,n);
    CVModule.predictorstep!(system,n);
    CVModule.correctorfluxes!(system,n);
    CVModule.correctorstep!(system,n);
    times.tfd += toq();
    tic();
    if hemoflag == "no"
        CVModule.applyendbcs!(system,n);
    elseif hemoflag == "yes"
        CVModule.applyendbcs!(system,n,hemoflag);
    end
    times.tc += toq();
    tic();
    CVModule.splitinvariants!(system,n);
    if hemoflag == "no"
        CVModule.solvesplits!(system,n);
    elseif hemoflag == "yes"
        CVModule.solvesplits!(system,n,hemoflag);
        # CVModule.applytourniquet!(system,n); # turn off to allow continual bleeding
    end
    times.ts += toq();
    tic();
    CVModule.arterialpressure!(system,n);
    CVModule.regulateall!(system,n);
    if coupleflag == "yes"
        CVModule.senddata(system,n,sender);
    end
    n+=1;
    times.tr += toq();
end
times.tt += toc();

CVModule.updatevolumes!(system,n);

if coupleflag == "yes"
    ZMQ.close(sender)
    ZMQ.close(ctx)
end

if saveflag == "yes"
    file = MAT.matopen("test3.mat", "w")
    write(file, "system", system)
    close(file)
end

return system, n, times

end
