importall CVModule
importall MAT

function main()

filename = "arterytree.csv";
# filename = "uc6.mat";
rstflag = "no"
hemoflag = "no"
saveflag = "yes"

system = CVModule.buildall(filename;numbeatstotal=1,restart=rstflag,injury=hemoflag);

# n = system.solverparams.nstart;
#
# tic();
# while system.solverparams.numbeats < system.solverparams.numbeatstotal
#     predictorfluxes!(system,n);
#     predictorstep!(system,n);
#     correctorfluxes!(system,n);
#     correctorstep!(system,n);
#     applyendbcs!(system,n);
#     splitinvariants!(system,n);
#     if hemoflag == "no"
#         solvesplits!(system,n);
#     elseif hemoflag == "yes"
#         solvesplits!(system,n,hemoflag);
#         applytourniquet!(system,n);
#     end
#     arterialpressure!(system,n);
#     regulateall!(system,n);
#     n+=1
# end
# toc()
#
# updatevolumes!(system,n);

if saveflag == "yes"
    file = MAT.matopen("test.mat", "w")
    write(file, "system", system)
    close(file)
end

return system

end
