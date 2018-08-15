importall CVModule

function main()

# options
loadfile = "arterytree.csv"; # default artery data file for new simulation
# loadfile = "run2.mat"; # filename for restarting simulation
rstflag = "no" # restarting from scratch or previous simulation
hemoflag = "no" # 10% hemorrhage from left femoral artery
saveflag = "yes" # save solution to .mat file
savefile = "run1.mat" # filename for saving (only used if saveflag == "yes")
coupleflag = "no" # coupling to 3D liver tissue model
timeflag = "yes" # solver timing

# build solution struct
system = CVModule.buildall(loadfile;numbeatstotal=1,restart=rstflag,injury=hemoflag);

# timers
times = CVModule.CVTimer();

# collect all IDs of terminal/non-terminal branches
terms = Int64[];
splits = Int64[];
for i = 1:length(system.branches.ID)
    if isempty(system.branches.children[i])
        push!(terms,i);
    else
        push!(splits,i);
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
    if mod(n,100) == 0
        println("Reached time step $n.")
    end
    if hemoflag == "no"
        CVModule.tvdrk3!(system,times,n,splits,terms);
        tic();
        CVModule.update0d!(system,n,terms);
    elseif hemoflag == "yes"
        CVModule.tvdrk3!(system,times,n,splits,terms,hemoflag);
        tic();
        CVModule.update0d!(system,n,terms,hemoflag);
    end
    times.tc += toq();
    tic();
    CVModule.arterialpressure!(system,n);
    CVModule.regulateall!(system,n,terms);
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
    file = MAT.matopen(savefile, "w")
    write(file, "system", system)
    close(file)
end

if timeflag == "yes"
    return system, n, times
elseif timeflag == "no"
    return system, n
end

end
