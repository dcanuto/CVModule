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

# filename = "arterytree.csv";
filename = "test1.mat";
rstflag = "yes"
hemoflag = "no"
saveflag = "yes"
coupleflag = "no"

# build solution struct
system = CVModule.buildall(filename;numbeatstotal=10,restart=rstflag,injury=hemoflag);

# timers
times = CVTimer();

# collect all IDs of terminal branches
terms = Int64[];
for i = 1:length(system.branches.ID)
    if isempty(system.branches.children[i])
        push!(terms,i);
    end
end

n = system.solverparams.nstart;

# coupling to 3D liver model
if coupleflag == "yes"
    ctx=ZMQ.Context();
    sender=ZMQ.Socket(ctx,ZMQ.REQ);
    ZMQ.connect(sender,"tcp://127.0.0.1:5555");
end

tic();
# solver loop
while system.solverparams.numbeats < system.solverparams.numbeatstotal
    tic();
    for i = 1:length(system.branches.ID)
        CVModule.predictorfluxes!(system.branches.Fp[i],system.branches.Q[i][:,n+1],
            system.branches.A[i][:,n+1],system.branches.beta[i][end],system.solverparams.rho,
            system.solverparams.acols,system.solverparams.qcols);
        CVModule.predictorstep!(system.branches.Aforward[i],system.branches.Qforward[i],
            system.branches.Abackward[i],system.branches.Qbackward[i],
            system.branches.A[i][:,n+1],system.branches.Q[i][:,n+1],
            system.solverparams.h,system.branches.k[i],system.branches.Fp[i],
            system.solverparams.colsint,system.solverparams.acolspre,system.solverparams.qcolspre,
            system.solverparams.diffusioncoeff,system.solverparams.mu,system.solverparams.rho);
        CVModule.correctorfluxes!(system.branches.Fbarforward[i],system.branches.Fbarbackward[i],
            system.branches.Qforward[i],system.branches.Aforward[i],system.branches.Qbackward[i],
            system.branches.Abackward[i],system.branches.beta[i][end],system.solverparams.rho,
            system.solverparams.acolscor,system.solverparams.qcolscor);
        ret1 = Array{Float64}(system.solverparams.JL-2);
        ret2 = Array{Float64}(system.solverparams.JL-2);
        CVModule.correctorstep!(ret1,ret2,
            system.branches.A[i][:,n+1],system.branches.Q[i][:,n+1],system.branches.Fbarforward[i],
            system.branches.Fbarbackward[i],system.branches.Qforward[i],system.branches.Aforward[i],
            system.branches.Qbackward[i],system.branches.Abackward[i],system.solverparams.colsint,
            system.solverparams.acolscor,system.solverparams.qcolscor,system.solverparams.h,system.branches.k[i],
            system.solverparams.diffusioncoeff,system.solverparams.mu,system.solverparams.rho);
        system.branches.A[i][system.solverparams.colsint,n+2] .= ret1;
        system.branches.Q[i][system.solverparams.colsint,n+2] .= ret2;
    end
    times.tfd += toq();
    tic();
    if hemoflag == "no"
        CVModule.applyendbcs!(system,n,terms);
    elseif hemoflag == "yes"
        @time CVModule.applyendbcs!(system,n,terms,hemoflag);
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
    file = MAT.matopen("test1.mat", "w")
    write(file, "system", system)
    close(file)
end

return system, n, times

end
