importall CVModule

function main()

rstflag = "no" # restart from prior solution, filename format must follow fnames below
hemoflag = "no" # 10% total blood vol. hemorrhage from left femoral artery
saveflag = "no" # save solution file to .mat struct
coupleflag = "no" # coupling with 3D organ model via ZMQ
assimflag = "yes" # patient data assimilation via EnKF to tune model state & params.
ptbflag = "yes" # generate ensemble via random perturbations, ONLY USE ONCE

ensemblesize = 2;
if rstflag == "no"
    fnames = ["arterytree.csv" for i=1:ensemblesize];
elseif rstflag == "yes"
    fnames = ["converge_$i.mat" for i=1:ensemblesize]
end

systems = pmap((a1)->CVModule.buildall(a1;numbeatstotal=1,restart=rstflag,
    injury=hemoflag,assim=assimflag),fnames);

if ptbflag == "yes"
    systems = pmap((a1)->CVModule.perturbics!(a1),systems);
end

n = [systems[1].solverparams.nstart for i=1:ensemblesize];

# if coupleflag == "yes"
#     ctx=ZMQ.Context();
#     sender=ZMQ.Socket(ctx,ZMQ.REQ);
#     ZMQ.connect(sender,"tcp://127.0.0.1:5555");
# end

# tic();
# while systems[1].solverparams.numbeats < systems[1].solverparams.numbeatstotal
#     soln = pmap((a1,a2)->CVModule.advancetime!(a1,a2;injury=hemoflag),systems,n);
#     systems = [soln[i][1] for i in 1:length(soln)];
#     n = [soln[i][2] for i in 1:length(soln)];
#     println("Current time: $(systems[1].t[n[1]+1])")
# end
# toc();

systems = pmap((a1,a2)->CVModule.updatevolumes!(a1,a2),systems,n);

# if coupleflag == "yes"
#     ZMQ.close(sender)
#     ZMQ.close(ctx)
# end

if saveflag == "yes"
    fnames = ["start_$i.mat" for i=1:ensemblesize];
    files = pmap((a1,a2)->MAT.matopen(a1,a2),fnames,["w" for i=1:ensemblesize]);
    pmap((a1,a2,a3)->write(a1,a2,a3),files,["system" for i=1:ensemblesize],systems);
    map((a1)->close(a1),files);
end

return systems, n

end
