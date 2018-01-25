importall CVModule

function main()

filename = "arterytree.csv";
# filename = "uc6.mat";
rstflag = "no"
hemoflag = "no"
saveflag = "yes"

system = CVModule.buildall(filename;numbeatstotal=1,restart=rstflag,injury=hemoflag);

n = system.solverparams.nstart;

tic();
while system.solverparams.numbeats < system.solverparams.numbeatstotal
    CVModule.predictorfluxes!(system,n);
    CVModule.predictorstep!(system,n);
    CVModule.correctorfluxes!(system,n);
    CVModule.correctorstep!(system,n);
    CVModule.applyendbcs!(system,n);
    CVModule.splitinvariants!(system,n);
    if hemoflag == "no"
        CVModule.solvesplits!(system,n);
    elseif hemoflag == "yes"
        CVModule.solvesplits!(system,n,hemoflag);
        CVModule.applytourniquet!(system,n);
    end
    CVModule.arterialpressure!(system,n);
    CVModule.regulateall!(system,n);
    n+=1
end
toc()

CVModule.updatevolumes!(system,n);

if saveflag == "yes"
    file = MAT.matopen("test.mat", "w")
    write(file, "system", system)
    close(file)
end

return system, n

end
