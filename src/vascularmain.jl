importall CVModule

function main()

filename = "arterytree.csv";
# filename = "test.mat";
rstflag = "no"
hemoflag = "no"
saveflag = "yes"
coupleflag = "no"
assimflag = "yes"

ensemblesize = 3;

fnames = [filename for i=1:ensemblesize];

systems = pmap((a1)->CVModule.buildall(a1;numbeatstotal=1,restart=rstflag,
    injury=hemoflag,assim=assimflag),fnames);

n = systems[1].solverparams.nstart;

# system = CVModule.buildall(filename;numbeatstotal=1,restart=rstflag,
#     injury=hemoflag,assim=assimflag);
#
# n = system.solverparams.nstart;

# if coupleflag == "yes"
#     ctx=ZMQ.Context();
#     sender=ZMQ.Socket(ctx,ZMQ.REQ);
#     ZMQ.connect(sender,"tcp://127.0.0.1:5555");
# end
#
# tic();
# while system.solverparams.numbeats < system.solverparams.numbeatstotal
#     CVModule.predictorfluxes!(system,n);
#     CVModule.predictorstep!(system,n);
#     CVModule.correctorfluxes!(system,n);
#     CVModule.correctorstep!(system,n);
#     CVModule.applyendbcs!(system,n);
#     CVModule.splitinvariants!(system,n);
#     if hemoflag == "no"
#         CVModule.solvesplits!(system,n);
#     elseif hemoflag == "yes"
#         CVModule.solvesplits!(system,n,hemoflag);
#         CVModule.applytourniquet!(system,n);
#     end
#     CVModule.arterialpressure!(system,n);
#     CVModule.regulateall!(system,n);
#     if coupleflag == "yes"
#         CVModule.senddata(system,n,sender);
#     end
#     if n > 0 && mod(n/system.pdata.nsamp,1) == 0
#         println(system.t[n+1])
#     end
#     n+=1;
# end
# toc();
#
# CVModule.updatevolumes!(system,n);
#
# if coupleflag == "yes"
#     ZMQ.close(sender)
#     ZMQ.close(ctx)
# end
#
# if saveflag == "yes"
#     file = MAT.matopen("save.mat", "w")
#     write(file, "system", system)
#     close(file)
# end

return systems, n

end
