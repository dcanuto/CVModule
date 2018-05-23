importall CVModule

function main()

filename = "arterytree.csv";
# filename = "test1.mat";
rstflag = "no"
hemoflag = "no"
saveflag = "yes"
coupleflag = "no"

system = CVModule.buildall(filename;numbeatstotal=1,restart=rstflag,injury=hemoflag);

n = system.solverparams.nstart;

if coupleflag == "yes"
    ctx=ZMQ.Context();
    sender=ZMQ.Socket(ctx,ZMQ.REQ);
    ZMQ.connect(sender,"tcp://127.0.0.1:5555");
end

tic();
while system.solverparams.numbeats < system.solverparams.numbeatstotal
    CVModule.predictorfluxes!(system,n);
    CVModule.predictorstep!(system,n);
    CVModule.correctorfluxes!(system,n);
    CVModule.correctorstep!(system,n);
    if hemoflag == "no"
        CVModule.applyendbcs!(system,n);
    elseif hemoflag == "yes"
        CVModule.applyendbcs!(system,n,hemoflag);
    end
    CVModule.splitinvariants!(system,n);
    if hemoflag == "no"
        CVModule.solvesplits!(system,n);
    elseif hemoflag == "yes"
        CVModule.solvesplits!(system,n,hemoflag);
        CVModule.applytourniquet!(system,n);
    end
    CVModule.arterialpressure!(system,n);
    CVModule.regulateall!(system,n);
    if coupleflag == "yes"
        CVModule.senddata(system,n,sender);
    end
    n+=1;
end
toc();

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

return system, n

end
