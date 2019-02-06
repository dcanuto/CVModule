importall CVModule

function main()

# options
rstflag = "yes" # restarting from scratch or previous simulation
hemoflag = "no" # 10% hemorrhage from left femoral artery
saveflag = "yes" # save solution to .mat file
coupleflag = "no" # coupling to 3D liver tissue model
timeflag = "yes" # solver timing
assimflag = "no" # patient data assimilation via EnKF
if rstflag == "no"
    sampleflag = "yes" # LH-OAT sampling (only use when startijng from scratch)
elseif rstflag == "yes"
    sampleflag = "no" # LH-OAT sampling (only use when starting from scratch)
end
colnum = 6;

# build solution struct or generate ensemble
if assimflag == "no"
    if rstflag == "no"
        loadfile = "arterytree.csv"; # default artery data file for new sim
    elseif rstflag == "yes"
        loadfile = "lhoat3_$colnum.mat"; # restart file
    end
    system = CVModule.buildall(loadfile;numbeatstotal=1,restart=rstflag,injury=hemoflag);
    savefile = "lhoat3_$colnum.mat" # filename for saving (only used if saveflag == "yes")
elseif assimflag == "yes"
    ensemblesize = 3;
    if rstflag == "no"
        loadfiles = ["arterytree.csv" for i=1:ensemblesize];
    elseif rstflag == "yes" loadfiles = ["test_1_$i.mat" for i=1:ensemblesize];
    end
    systems = pmap((a1)->CVModule.buildall(a1;numbeatstotal=1,restart=rstflag,
        injury=hemoflag),loadfiles);
    savefiles = ["test_1_$i.mat" for i=1:ensemblesize];
end

# ensemble distributions, allocators, and scalings
if assimflag == "yes"
    nparams = 1;
    nstates = 1;
    nmeas = 1;
    errors = CVModule.Errors(nparams);

    X = [zeros(nstates) for i=1:ensemblesize];
    θ = [zeros(nparams) for i=1:ensemblesize];
    Y = [zeros(nmeas) for i=1:ensemblesize];

    xhat = zeros(nstates);
    θhat = zeros(nparams);
    yhat = zeros(nmeas);

    yi = [zeros(nmeas) for i=1:ensemblesize];

    θs = ones(nparams);
    for i=1:ensemblesize
        θs[1] += 1.;
    end
    θs /= ensemblesize;

    p = 0.5; # RTPS relaxation amount, ∃ [0.5,0.95]
    c = zeros(ensemblesize);
    σxb = zeros(nstates);
    σxa = zeros(nstates);
    σtb = zeros(nparams);
    σta = ones(nparams);

    tout = Float64[];
    xout = Vector{Float64}[];
    xoutv = Vector{Float64}[];
    yout = Vector{Float64}[];
    youtv = Vector{Float64}[];
    θout = Vector{Float64}[];
    θoutv = Vector{Float64}[];
    Pθout = Vector{Float64}[];
    Pθoutv = Vector{Float64}[];
    Pxout = Vector{Float64}[];
    Pxoutv = Vector{Float64}[];
    lbx = Vector{Float64}[];
    ubx = Vector{Float64}[];
    lbxv = Vector{Float64}[];
    ubxv = Vector{Float64}[];
end

# timers
if assimflag == "no"
    times = CVModule.CVTimer();
elseif assimflag == "yes"
    times = [CVModule.CVTimer() for i=1:ensemblesize];
end

# collect all IDs of terminal/non-terminal branches
terms = Int64[];
splits = Int64[];
if assimflag == "no"
    for i = 1:length(system.branches.ID)
        if isempty(system.branches.children[i])
            push!(terms,i);
        else
            push!(splits,i);
        end
    end
elseif assimflag == "yes"
    for i = 1:length(systems[1].branches.ID)
        if isempty(systems[1].branches.children[i])
            push!(terms,i);
        else
            push!(splits,i);
        end
    end
    term_itr = [terms for i=1:ensemblesize];
    split_itr = [splits for i=1:ensemblesize];
end

# time step counter
if assimflag == "no"
    n = system.solverparams.nstart;
elseif assimflag == "yes"
    n = [systems[1].solverparams.nstart for i=1:ensemblesize];
end

# insert LH-OAT sample into model parameters
if sampleflag == "yes"
    rstflag == "no" ||
        throw(ArgumentError("rstflag must be set to \"no\" to insert LH-OAT sample"))
    pointfile = "lhoat_3.mat"
    vname = "lhoat_3"
    pin = MAT.matread(pointfile);
    params = pin[vname][:,colnum];
    CVModule.sampletoparams!(system,params,rstflag,assimflag,sampleflag)
end

# coupling to 3D liver model (only for non-ensemble simulation)
if coupleflag == "yes" && assimflag == "no"
    ctx=ZMQ.Context();
    sender=ZMQ.Socket(ctx,ZMQ.REQ);
    ZMQ.connect(sender,"tcp://127.0.0.1:5555");
end

# solver loop
tic();
if assimflag == "no"
    tic();
    while system.solverparams.numbeats < system.solverparams.numbeatstotal
        # if mod(n,100) == 0
        #     println("Reached time step $n.")
        # end
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
elseif assimflag == "yes"
    nsamp = [1000 for i=1:ensemblesize];
    while systems[1].solverparams.numbeats < systems[1].solverparams.numbeatstotal
        # normal time integration between measurements
        rflag = "predict";
        soln = pmap((a1,a2,a3,a4,a5,a6)->CVModule.advancetime!(a1,a2,a3,a4,a5,a6;
            injury=hemoflag,runtype=rflag),systems,times,n,nsamp,term_itr,split_itr);
        systems = [soln[i][1] for i=1:ensemblesize];
        n = [soln[i][2] for i=1:ensemblesize];
        times = [soln[i][3] for i=1:ensemblesize];
        # data assimilation
        if assimflag == "yes" && mod((n[1]-1)/nsamp[1],1) == 0
            println("Current time step: $(n[1]-1)")
            println("Current time: $(systems[1].t[n[1]])")

            # # non-dimensional measurements (TODO)
            # y = zeros(nmeas);
            #
            # # generate non-dimensional measurement replicates
            # for i = 1:ensemblesize
            #     for j = 1:nmeas
            #         yi[i][j] = y[j] + rand(Distributions.Normal(0,error.odev[j]/error.oscales[j]));
            #     end
            # end
            #
            # # forecast parameters and means
            # for i = 1:ensemblesize
            #     θ[i][1] = 0.;
            # end
            #
            # # forecast parameters (dimensional)
            # θ = [CVModule.paramwalk!(error,θs,θ[i]) for i=1:ensemblesize];
            #
            # # non-dimensionalize parameters
            # for i = 1:ensemblesize
            #     θ[i] = θ[i]./θs;
            # end
            #
            # # RTPS prior std. dev.
            # for i = 1:nparams
            #     for j = 1:ensemblesize
            #         c[j] = θ[j][i];
            #     end
            #     σtb[i] = std(c;corrected=true);
            # end
            #
            # # parameters back into model
            # for i = 1:ensemblesize
            # end

            # calculate elastance scaling (TODO)

            # forecast state w/ forecast parameters (single time step)
            rflag = "assim";
            for i = 1:ensemblesize
                systems[i],n[i],times[i] = CVModule.advancetime!(systems[i],
                    times[i],n[i],nsamp[i],term_itr[i],split_itr[i];
                    injury=hemoflag,runtype=rflag)
            end
            # soln = pmap((a1,a2,a3,a4,a5,a6)->CVModule.advancetime!(a1,a2,a3,a4,a5,a6;
            #     injury=hemoflag,runtype=rflag),systems,times,n,nsamp,term_itr,split_itr);
            # systems = [soln[i][1] for i=1:ensemblesize];
            # n = [soln[i][2] for i=1:ensemblesize];
            # times = [soln[i][3] for i=1:ensemblesize];

            # # vector of forecast measurements
            # Y = [zeros(nmeas) for i=1:ensemblesize];
            #
            # # forecast mean params, measurements
            # yhat = mean(Y);
            # θhat = mean(θ);
            #
            # # covariance/cross-covariance
            # Pty = zeros(nparams,nmeas);
            # Pyy = zeros(nmeas,nmeas);
            # for i=1:ensemblesize
            #     Pty += *((θ[i] .- θhat),(Y[i] .- yhat)');
            #     Pyy += *((Y[i] .- yhat),(Y[i] .- yhat)');
            # end
            # Pty ./= ensemblesize;
            # Pyy ./= ensemblesize;
            #
            # # add noise to meas. covariance (allows invertibility)
            # Pyy += diagm((error.odev./error.oscale)^2,0);
            #
            # # parameter Kalman gain
            # K = Pty*inv(Pyy);
            #
            # # parameter analysis step
            # for i = 1:ensemblesize
            #     θ[i][:] += K*(yi[i] .- Y[i]);
            # end
            #
            # # RTPS parameter covariance inflation
            # for i = 1:nparams
            #     for j = 1:ensemblesize
            #         c[j] = θ[j][i];
            #     end
            #     σta = std(c;corrected=true);
            # end
            # θhat = mean(θ);
            # for i = 1:ensemblesize
            #     θ[i] .= ϴ[i] .+ p.*((σtb.-σta)./σta).*(θ[i].-θhat);
            # end

            # analysis parameters back into ensemble
            for i = 1:ensemblesize
                Rdefault = [0.3,0.21,0.003,0.01]*CVModule.mmHgToPa/CVModule.cm3Tom3;
                Cdefault = [0.01,1.64,1.81,13.24,73.88]*CVModule.cm3Tom3/CVModule.mmHgToPa;
                CVModule.assignterminals!(systems[i],Rdefault,Cdefault,Dict("a"=>0),"no","yes")
            end

            # # recalculate parameter distribution
            # θhat = mean(θ);
            # append!(θout,[θhat])
            #
            # Ptt = zeros(nparams);
            # for i = 1:ensemblesize
            #     Ptt += (θ[i] .- θhat).^2;
            # end
            # Ptt ./= ensemblesize;
            # append!(Pθout,[Ptt])

            # calculate elastance scaling (TODO)

            # corrected forecast w/ analysis parameters
            rflag = "assim";
            for i = 1:ensemblesize
                systems[i],n[i],times[i] = CVModule.advancetime!(systems[i],
                    times[i],n[i],nsamp[i],term_itr[i],split_itr[i];
                    injury=hemoflag,runtype=rflag)
            end
            # soln = pmap((a1,a2,a3,a4,a5,a6)->CVModule.advancetime!(a1,a2,a3,a4,a5,a6;
            #     injury=hemoflag,runtype=rflag),systems,times,n,nsamp,term_itr,split_itr);
            # systems = [soln[i][1] for i=1:ensemblesize];
            # n = [soln[i][2] for i=1:ensemblesize];
            # times = [soln[i][3] for i=1:ensemblesize];

        #     # forecast state, measurements
        #     X = [zeros(nstates) for i=1:ensemblesize];
        #     Y = [zeros(nmeas) for i=1:ensemblesize];
        #
        #     # mean state, measurement
        #     xhat = mean(X);
        #     yhat = mean(Y);
        #
        #     # output means
        #     append!(xout,[xhat])
        #     append!(yout,[yhat])
        #
        #     # state variance, 2-sigma quantiles (for post-processing)
        #     Pxx = zeros(length(xhat));
        #     for i = 1:ensemblesize
        #         Pxx += (X[i] .- xhat).^2;
        #     end
        #     Pxx ./= ensemblesize;
        #     append!(Pxout,[Pxx])
        #
        #     for i = 1:nstates
        #         xq = Float64[];
        #         for j = 1:ensemblesize
        #             push!(xq,X[j][i])
        #         end
        #         q = quantile(xq,[0.025,0.975]);
        #         if i == 1
        #             append!(lbx,[ones(1)*q[1]])
        #             append!(ubx,[ones(1)*q[2]])
        #         else
        #             push!(lbx[end],q[1])
        #             push!(ubx[end],q[2])
        #         end
        #     end
        end

        # # output measurement time
        # push!(tout,systems[1].t[n[1]])
    end
end
toc();

if assimflag == "no"
    CVModule.updatevolumes!(system,n);
elseif assimflag == "yes"
    for i = 1:ensemblesize
        systems[i] = CVModule.updatevolumes!(systems[i],n[i])
    end
    # systems = pmap((a1,a2)->CVModule.updatevolumes!(a1,a2),systems,n);
end

# # reshape ensemble averages into vectors of individual time series
# if assimflag == "yes"
#     for i in 1:length(xout[1]) # i indexes state variables
#         xo = Float64[];
#         Pxo = Float64[];
#         lb = Float64[];
#         ub = Float64[];
#         for j in 1:length(xout) # j indexes time steps
#             push!(xo,xout[j][i])
#             push!(Pxo,Pxout[j][i])
#             push!(lb,lbx[j][i])
#             push!(ub,ubx[j][i])
#         end
#         append!(xoutv,[xo])
#         append!(Pxoutv,[Pxo])
#         append!(lbxv,[lb])
#         append!(ubxv,[ub])
#     end
#
#     for i in 1:length(θout[1]) # i indexes state variables
#         θo = Float64[];
#         Pθo = Float64[];
#         for j in 1:length(θout) # j indexes time steps
#             push!(θo,θout[j][i])
#             push!(Pθo,Pθout[j][i])
#         end
#         append!(θoutv,[θo])
#         append!(Pθoutv,[Pθo])
#     end
#
#     for i in 1:length(yout[1])
#         yo = Float64[];
#         for j in 1:length(yout)
#             push!(yo,yout[j][i])
#         end
#         append!(youtv,[yo])
#     end
# end

if coupleflag == "yes" && assimflag == "no"
    ZMQ.close(sender)
    ZMQ.close(ctx)
end

if saveflag == "yes"
    if assimflag == "no"
        file = MAT.matopen(savefile, "w")
        write(file,"system",system)
        close(file)
    elseif assimflag == "yes"
        for i=1:ensemblesize
            file = MAT.matopen(savefiles[i],"w");
            write(file,"system",systems[i])
            close(file)
        end
    end
end

if assimflag == "no"
    if timeflag == "yes"
        return system, n, times
    elseif timeflag == "no"
        return system, n
    end
elseif assimflag == "yes"
    return n,times,term_itr,split_itr
end

end
