importall CVModule

function main()

filename = "arterytree.csv";
# filename = "test.mat";
rstflag = "no"
hemoflag = "no"
saveflag = "yes"
coupleflag = "no"
assimflag = "yes"

ensemblesize = 2;

fnames = [filename for i=1:ensemblesize];

systems = pmap((a1)->CVModule.buildall(a1;numbeatstotal=1,restart=rstflag,
    injury=hemoflag,assim=assimflag),fnames);

n = [systems[1].solverparams.nstart for i=1:ensemblesize];

# if coupleflag == "yes"
#     ctx=ZMQ.Context();
#     sender=ZMQ.Socket(ctx,ZMQ.REQ);
#     ZMQ.connect(sender,"tcp://127.0.0.1:5555");
# end

tic();
systems = pmap((a1,a2)->CVModule.advancetime!(a1,a2;injury=hemoflag),systems,n);
toc();

systems = pmap((a1,a2)->CVModule.updatevolumes!(a1,a2),systems,n);

# if coupleflag == "yes"
#     ZMQ.close(sender)
#     ZMQ.close(ctx)
# end
#
if saveflag == "yes"
    fnames = ["save_$i.mat" for i=1:ensemblesize];
    files = pmap((a1,a2)->MAT.matopen(a1,a2),fnames,["w" for i=1:ensemblesize]);
    pmap((a1,a2,a3)->write(a1,a2,a3),files,["system" for i=1:ensemblesize],systems);
    pmap((a1)->close(a1),files);
end

return systems, n

end
