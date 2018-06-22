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

# options
filename = "arterytree.csv"; # default artery data file for new simulation
# filename = "test1.mat"; # filename for restarting simulation
rstflag = "no" # restarting from scratch or previous simulation
hemoflag = "no" # 10% hemorrhage from left femoral artery
saveflag = "yes" # save solution to .mat file
coupleflag = "no" # coupling to 3D liver tissue model
timeflag = "yes" # solver timing

# build solution struct
system = CVModule.buildall(filename;numbeatstotal=1,restart=rstflag,injury=hemoflag);

# timers
times = CVTimer();

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
ret1 = Array{Float64}(system.solverparams.JL-2); # placeholder arrays for interior 1d points
ret2 = Array{Float64}(system.solverparams.JL-2);
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
        CVModule.applyendbcs!(system,n,terms,hemoflag);
    end
    times.tc += toq();
    tic();
    for i = 1:length(splits)
        children = system.branches.children[splits[i]];
        if length(children) == 1
            f = CVModule.fsingle;
            J = CVModule.Jsingle;
        elseif length(children) == 2
            f = CVModule.fdouble;
            J = CVModule.Jdouble;
        elseif length(children) == 3
            f = CVModule.ftriple;
            J = CVModule.Jtriple;
        elseif length(children) == 4
            f = CVModule.fquad;
            J = CVModule.Jquad;
        end
        W = zeros(length(children)+1);
        beta = zeros(length(children)+1);
        A0 = zeros(length(children)+1);
        c0 = zeros(length(children)+1);
        k = zeros(length(children)+1);
        Q = zeros(length(children)+1);
        Qo = zeros(length(children)+1);
        A = zeros(length(children)+1);
        Ao = zeros(length(children)+1);
        beta[1] = system.branches.beta[splits[i]][end];
        A0[1] = system.branches.A0[splits[i]][end];
        c0[1] = system.branches.c0[splits[i]][end];
        k[1] = system.branches.k[splits[i]][end];
        Q[1] = system.branches.Q[splits[i]][system.solverparams.JL,n+1];
        Qo[1] = system.branches.Q[splits[i]][system.solverparams.JL-1,n+1];
        A[1] = system.branches.A[splits[i]][system.solverparams.JL,n+1];
        Ao[1] = system.branches.A[splits[i]][system.solverparams.JL-1,n+1];
        for j = 1:length(children)
            beta[j+1] = system.branches.beta[children[j]][end];
            A0[j+1] = system.branches.A0[children[j]][end];
            c0[j+1] = system.branches.c0[children[j]][end];
            k[j+1] = system.branches.k[children[j]][end];
            Q[j+1] = system.branches.Q[children[j]][1,n+1];
            Qo[j+1] = system.branches.Q[children[j]][2,n+1];
            A[j+1] = system.branches.A[children[j]][1,n+1];
            Ao[j+1] = system.branches.A[children[j]][2,n+1];
        end
        CVModule.splitinvariants!(W,Q,A,Qo,Ao,beta,c0,system.solverparams.rho,
            k,system.solverparams.h,system.solverparams.diffusioncoeff,
            system.solverparams.mu,children);
        system.branches.W1[splits[i]] = W[1];
        system.branches.W2[children] .= W[2:end];
        Qnew = zeros(length(children)+1);
        Anew = zeros(length(children)+1);
        iters = Int64[0];
        if hemoflag == "no"
            CVModule.solvesplits!(iters,Qnew,Anew,children,Q,A,W,beta,A0,c0,
                system.solverparams.rho,f,J,system.solverparams.maxiter,
                system.solverparams.maxval,system.solverparams.epsJ,
                system.solverparams.epsN)
        elseif hemoflag == "yes"
            curloss = [0.];
            CVModule.solvesplits!(iters,Qnew,Anew,children,Q,A,W,beta,A0,c0,
                system.solverparams.rho,f,J,system.solverparams.maxiter,
                system.solverparams.maxval,system.solverparams.epsJ,
                system.solverparams.epsN,hemoflag,system.hemo.injured[splits[i]],
                system.solverparams.h,system.hemo.Ph,curloss);
            system.hemo.Vloss += curloss[1];
            if system.hemo.injured[splits[i]] == true
                system.branches.W2[splits[i]] = W[1];
                system.branches.W1[children] .= W[2:end];
            end
            CVModule.applytourniquet!(system,n); # turn off to allow continual bleeding
        end
        system.solverparams.totaliter += iters[1];
        system.branches.Q[splits[i]][system.solverparams.JL,n+2] = Qnew[1];
        system.branches.A[splits[i]][system.solverparams.JL,n+2] = Anew[1];
        for j = 1:length(children)
            system.branches.Q[children[j]][1,n+2] = Qnew[j+1];
            system.branches.A[children[j]][1,n+2] = Anew[j+1];
        end
    end
    times.ts += toq();
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
    file = MAT.matopen("savefile.mat", "w")
    write(file, "system", system)
    close(file)
end

if timeflag == "yes"
    return system, n, times
elseif timeflag == "no"
    return system, n
end

end
