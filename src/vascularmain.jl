importall CVModule

function main()

rstflag = "no" # restart from prior solution, filename format must follow fnames below
hemoflag = "no" # 10% total blood vol. hemorrhage from left femoral artery
saveflag = "yes" # save solution file to .mat struct
coupleflag = "no" # coupling with 3D organ model via ZMQ
assimflag = "yes" # patient data assimilation via EnKF to tune model state & params.
ptbflag = "no" # generate ensemble via random perturbations, ONLY USE ONCE

ensemblesize = 6;
if rstflag == "no"
    fnames = ["arterytree.csv" for i=1:ensemblesize];
elseif rstflag == "yes"
    fnames = ["p30_1_$i.mat" for i=1:ensemblesize];
end

systems = pmap((a1)->CVModule.buildall(a1;numbeatstotal=1,restart=rstflag,
    injury=hemoflag,assim=assimflag),fnames);

# if ptbflag == "yes"
#     for i = 1:length(systems) # using 1st system as reference for measurements
#         CVModule.perturbics!(systems[i])
#         # println("β for ensemble member $(i): $(systems[i].branches.beta[61][end])")
#         # println("A0 for ensemble member $(i): $(systems[i].branches.A0[61][end])")
#     end
# end

n = [systems[1].solverparams.nstart for i=1:ensemblesize];

# if coupleflag == "yes"
#     ctx=ZMQ.Context();
#     sender=ZMQ.Socket(ctx,ZMQ.REQ);
#     ZMQ.connect(sender,"tcp://127.0.0.1:5555");
# end

if assimflag == "yes"
    numarteries = 4;
    # allocators for ensemble augmented state, measurements
    # X = [zeros(numarteries*(systems[1].solverparams.JL + length(systems[1].error.pdev)-1)+3) for i in (1:length(systems))];
    X = [zeros(numarteries*(systems[1].solverparams.JL)+2) for i in (1:length(systems))];
    θ = [zeros(length(systems[1].error.pdev)) for i in (1:length(systems))];

    # Y = [zeros(systems[1].solverparams.JL) for i in (1:length(systems)-1)];
    Y = [zeros(1) for i in (1:length(systems))];

    # allocators for state, forecast measurement mean
    # xhat = zeros(numarteries*(systems[1].solverparams.JL + length(systems[1].error.pdev)-1)+3);
    xhat = zeros(numarteries*(systems[1].solverparams.JL)+2);
    θhat = zeros(length(systems[1].error.pdev));
    yhat = zeros(1);

    # allocators for measurement replicates
    yi = zeros(length(systems));

    # index variables for state vector
    startindex = 0;
    endindex = 0;

    # output variables
    println("Number of steps: $(systems[1].solverparams.numsteps)")
    println("Steps between samples: $(systems[1].pdata.nsamp)")
    βhat = zeros(Int(cld(systems[1].solverparams.numsteps,systems[1].pdata.nsamp)),numarteries);
    A0hat = zeros(Int(cld(systems[1].solverparams.numsteps,systems[1].pdata.nsamp)),numarteries);
    Phat = zeros(Int(cld(systems[1].solverparams.numsteps,systems[1].pdata.nsamp)));
    τ1hat = zeros(Int(cld(systems[1].solverparams.numsteps,systems[1].pdata.nsamp)));
    numassims = 0;
    println("Number of assimilations: $(Int(cld(systems[1].solverparams.numsteps,systems[1].pdata.nsamp)))")

    # allocators for calculating normalized RMSE ratio
    r1dot = 0;
    r2dot = zeros(ensemblesize);

    # load patient arterial pressure, scale to model HR
    filename = "pdata.mat"
    for i = 1:length(systems)
        scalepdata!(systems[i],filename);
    end

    # non-dimensional scalings
    As = systems[1].branches.A[61][1,6];
    Vs = systems[1].heart.lv.V[1];
    Ps = maximum(systems[1].pdata.P);
    θs = zeros(length(θhat));

    # RTPS allocators
    p = 0.95; # relaxation amount, ∃ [0.5, 0.95]
    c = zeros(length(systems));
    σb = zeros(length(xhat));
    σa = zeros(length(xhat));
    σtb = zeros(length(θhat));
    σta = zeros(length(θhat));

    # discretized times for determining ventricular elastance scaling
    t = linspace(0,systems[1].heart.activation.th[end],10000);
end

tic();
while systems[1].solverparams.numbeats < systems[1].solverparams.numbeatstotal
    if n[1] == 0 && rstflag == "no"
        Pd = 50; # estimated diastolic pressure, mmHg
        Vdl = systems[1].heart.lv.V[1]/CVModule.cm3Tom3;
        Vdr = systems[1].heart.lv.V[1]/CVModule.cm3Tom3;
        for j = 1:length(systems)
            for i = 61:64
                # randomize initial areas to avoid σ = 0 for arterial area
                Ps = rand(Distributions.Normal(Pd,1));
                systems[j].branches.A[i][1,:] = (Ps*CVModule.mmHgToPa/systems[j].branches.beta[i][end] +
                    sqrt(systems[j].branches.A0[i][end]))^2;
            end
            Vsl = rand(Distributions.Normal(Vdl,2.5));
            Vsr = rand(Distributions.Normal(Vdr,2.5));
            systems[j].heart.lv.V[1] = Vsl*CVModule.cm3Tom3;
            systems[j].heart.rv.V[1] = Vsr*CVModule.cm3Tom3;
        end
    end
    rflag = "predict";
    soln = pmap((a1,a2)->CVModule.advancetime!(a1,a2;injury=hemoflag,runtype=rflag),systems,n);
    systems = [soln[i][1] for i in 1:length(soln)];
    n = [soln[i][2] for i in 1:length(soln)];
    if assimflag == "yes" && mod((n[1]-1)/systems[1].pdata.nsamp,1) == 0
        rflag = "assim";

        println("Current time step: $(n[1])")
        println("Current time: $(systems[1].t[n[1]+1])")
        println("Measured pressure: $(systems[1].pdata.P[numassims+1])")

        # measurement
        y = systems[1].pdata.P[numassims+1];
        y /= Ps;

        # generate measurement replicates
        for i = 1:length(systems)
            yi[i] = y + rand(Distributions.Normal(0,systems[i].error.odev[1]/Ps));
        end

        # forecast parameters
        for i = 1:length(systems)
            for j = 1:numarteries
                θ[i][j] = systems[i].branches.beta[60+j][end];
            end
            θ[i][end] = systems[i].heart.activation.tau1;
        end

        # parameter means
        for j = 61:64
            θs[j-60] = 0;
            for i = 1:length(soln)
                θs[j-60] += systems[i].branches.beta[j][end];
            end
            θs[j-60] /= (ensemblesize);
            for i = 1:length(soln)
                systems[i].error.pbar[j-60] = θs[j-60];
            end
        end
        θs[end] = 0;
        for i = 1:length(soln)
            θs[end] += systems[i].heart.activation.tau1;
        end
        θs[end] /= (ensemblesize);
        for i = 1:length(soln)
            systems[i].error.pbar[end] = θs[end];
        end

        # forecast parameters
        θ = [CVModule.paramwalk!(systems[i],θ[i]) for i in 1:length(systems)];
        for i = 1:length(systems)
            # println("Forecast parameters, member $i: $(θ[i])")
            θ[i] = θ[i]./θs;
        end

        # RTPS prior standard deviation
        for i in 1:length(θ[1])
            for j in 1:length(systems)
                c[j] = θ[j][i];
            end
            σtb[i] = std(c;corrected=true);
        end

        # parameters back into model
        for i = 1:length(soln)
            for j = 1:numarteries
                systems[i].branches.beta[60+j][end] = θ[i][j]*θs[j];
            end
            systems[i].heart.activation.tau1 = θ[i][end]*θs[end];
            g1 = (t/systems[i].heart.activation.tau1).^systems[i].heart.activation.m1;
            g2 = (t/systems[i].heart.activation.tau2).^systems[i].heart.activation.m2;
            systems[i].heart.activation.k[end] = maximum((g1./(1+g1)).*(1./(1+g2)))^-1;
        end

        # forecast state
        soln = map((a1,a2)->CVModule.advancetime!(a1,a2;injury=hemoflag,runtype=rflag),systems,n);
        systems = [soln[i][1] for i in 1:length(soln)];
        X = [systems[i].branches.A[61][n[i]+1,:]/As for i in 1:length(soln)];
        for j = 62:64
            for i = 1:length(soln)
                append!(X[i],systems[i].branches.A[j][n[i]+1,:]/As)
            end
        end
        for i = 1:length(soln)
            append!(X[i],[systems[i].heart.lv.V[n[i]+1],systems[i].heart.rv.V[n[i]+1]]/Vs)
        end
        # println("Size of X, first forecast: $(size(X))")

        # forecast measurements
        Y = [systems[i].branches.P[61][n[i]+1,6]/Ps for i in 1:length(soln)];
        # println("Type of Y: $(typeof(Y))")
        # println("Size of Y: $(size(Y))")

        # forecast mean state, parameters, measurement
        xhat = mean(X);
        yhat = mean(Y);
        θhat = mean(θ);
        # println("Size of ̂x: $(size(xhat))")
        println("Normalized ̂x, first forecast: $xhat")
        println("Normalized ̂y, first forecast: $yhat")
        println("Normalized ̂θ, first forecast: $θhat")

        # forecast params./meas. cross covariance, meas. covariance
        Pty = zeros(length(θhat),length(yhat))
        Pyy = zeros(length(yhat),length(yhat))
        for i = 1:length(soln)
            Pty += *((θ[i] - θhat),(Y[i] - yhat)');
            Pyy += *((Y[i] - yhat),(Y[i] - yhat)');
        end
        Pty /= (ensemblesize);
        Pyy /= (ensemblesize);
        # add meas. noise to meas. covariance (allows invertibility)
        Pyy += diagm((systems[1].error.odev[1]/Ps)^2*ones(length(yhat)),0);

        # parameter Kalman gain
        K = Pty*inv(Pyy);
        # println("Kθ: $K")
        # println("Type of K: $(typeof(K))")
        # println("Size of K: $(size(K))")

        # parameter analysis step
        for i = 1:length(soln)
            # println("Type of correction: $(typeof(K*(yi - Y[i-1])))")
            # println("Size of correction: $(size(K*(yi - Y[i-1])))")
            θ[i][:] += K*(yi[i] - Y[i]);
            # println("Normalized analysis τ1, ensemble member $i: $(θ[i][end])")
        end
        # println("Normalized analysis ̂θ before RTPS inflation: $(mean(θ))")

        # RTPS parameter covariance inflation
        for i in 1:length(θ[1])
            for j in 1:length(systems)
                c[j] = θ[j][i];
            end
            σta[i] = std(c;corrected=true);
        end
        θhat = mean(θ);
        for i = 1:length(soln)
            θ[i][:] = θ[i][:] + p*((σtb-σta)./σta).*(θ[i][:]-θhat);
        end
        # println("Normalized analysis ̂θ after RTPS inflation: $(mean(θ))")
        # println("σθ difference = 0: $([(σtb[i]-σta[i]) == 0 for i in 1:length(σtb)])")

        # analysis parameters back into ensemble members
        for i = 1:length(soln)
            # systems[i].error.pbar[:] = mean(θ).*θs;
            for j = 1:numarteries
                systems[i].branches.beta[60+j][end] = θ[i][end-(numarteries-j)-1]*θs[j];
                if systems[i].branches.beta[60+j][end] <= systems[i].error.lb[j]
                    println("Warning: analysis β non-positive for artery $(60+j), member $i.
                        Setting to lower bound of $(systems[i].error.lb[j]) Pa/m.")
                    systems[i].branches.beta[60+j][end] = systems[i].error.lb[j];
                end
            end
            systems[i].heart.activation.tau1 = θ[i][end]*θs[end];
            if systems[i].heart.activation.tau1 <= systems[i].error.lb[end]
                println("Warning: analysis τ1 non-positive for member $i.
                    Setting to lower bound of $(systems[i].error.lb[end]) s.")
                systems[i].heart.activation.tau1 = systems[i].error.lb[end];
            end
            g1 = (t/systems[i].heart.activation.tau1).^systems[i].heart.activation.m1;
            g2 = (t/systems[i].heart.activation.tau2).^systems[i].heart.activation.m2;
            systems[i].heart.activation.k[end] = maximum((g1./(1+g1)).*(1./(1+g2)))^-1;
            # println("Analysis β, ensemble member $(i): $(systems[i].branches.beta[61][end])")
        end
        for i = 1:length(soln)
            for j = 1:numarteries
                βhat[numassims+1,j] += systems[i].branches.beta[60+j][end];
            end
            τ1hat[numassims+1] += systems[i].heart.activation.tau1;
        end
        βhat[numassims+1,:] /= ensemblesize;
        τ1hat[numassims+1] /= ensemblesize;
        for i = 1:length(soln)
            systems[i].error.pbar[1:end-1] = βhat[numassims+1,:];
            systems[i].error.pbar[end] = τ1hat[numassims+1];
        end
        println("Analysis mean β: $(βhat[numassims+1,:])")
        println("Analysis mean τ1: $(τ1hat[numassims+1])")

        # corrected forecast w/ analysis parameters
        soln = map((a1,a2)->CVModule.advancetime!(a1,a2;injury=hemoflag,runtype=rflag),systems,n);
        systems = [soln[i][1] for i in 1:length(soln)];
        X = [systems[i].branches.A[61][n[i]+1,:]/As for i in 1:length(soln)];
        for j = 62:64
            for i = 1:length(soln)
                append!(X[i],systems[i].branches.A[j][n[i]+1,:]/As)
            end
        end
        for i = 1:length(soln)
            append!(X[i],[systems[i].heart.lv.V[n[i]+1],systems[i].heart.rv.V[n[i]+1]]/Vs)
        end
        # println("Size of X, second forecast: $(size(X))")

        # prior standard deviation
        for i in 1:length(X[1])
            for j in 1:length(systems)
                c[j] = X[j][i];
            end
            σb[i] = std(c;corrected=true);
        end

        # second forecast mean state, measurement
        xhat = mean(X);
        yhat = mean(Y);
        # println("Size of ̂x: $(size(xhat))")
        # Phat[numassims+1] = yhat[end]*Ps;
        println("Normalized ̂x, second forecast: $xhat")
        println("Normalized ̂y, second forecast: $yhat")

        # forecast state/meas. cross covariance, meas. covariance
        Pxy = zeros(length(xhat),length(yhat))
        Pyy = zeros(length(yhat),length(yhat))
        for i = 1:length(soln)
            Pxy += *((X[i] - xhat),(Y[i] - yhat)');
            Pyy += *((Y[i] - yhat),(Y[i] - yhat)');
        end
        Pxy /= (ensemblesize);
        Pyy /= (ensemblesize);
        # add meas. noise to meas. covariance (allows invertibility)
        Pyy += diagm((systems[1].error.odev[1]/Ps)^2*ones(length(yhat)),0);

        # Kalman gain
        K = Pxy*inv(Pyy);
        # println("Kx: $K")
        # println("Type of K: $(typeof(K))")
        # println("Size of K: $(size(K))")

        # analysis step, NRR tracking
        for i = 1:length(soln)
            # println("Type of correction: $(typeof(K*(yi - Y[i-1])))")
            # println("Size of correction: $(size(K*(yi - Y[i-1])))")
            X[i][:] += K*(yi[i] - Y[i]);
            r2dot[i] += dot((Y[i]-yi),(Y[i]-yi));
        end
        r1dot += sqrt.(dot((yhat-y),(yhat-y)));

        # RTPS multiplicative covariance inflation
        for i in 1:length(X[1])
            for j in 1:length(systems)
                c[j] = X[j][i];
            end
            σa[i] = std(c;corrected=true);
        end
        xhat = mean(X);
        # println("Prior standard deviation: $(σb)")
        # println("Posterior standard deviation: $(σa)")
        # println("Standard deviation difference: $(σb-σa)")
        # println("Normalized difference: $((σb-σa)./σa)")
        # println("σx difference = 0: $([(σb[i]-σa[i]) == 0 for i in 1:length(σb)])")
        for i = 1:length(soln)
            X[i][:] = X[i][:] + p*((σb-σa)./σa).*(X[i][:]-xhat);
            # println("Analysis X, member $i: $(X[i])")
        end

        # analysis back into ensemble members
        for i = 1:length(soln)
            for j = 1:numarteries
                startindex = systems[i].solverparams.JL*(j-1)+1;
                endindex = systems[i].solverparams.JL*j;
                systems[i].branches.A[60+j][n[i]+1,:] = X[i][startindex:endindex]*As;
            end
            systems[i].heart.lv.V[n[i]+1] = X[i][endindex+1]*Vs;
            systems[i].heart.rv.V[n[i]+1] = X[i][endindex+2]*Vs;
        end

        # recompute ventricular pressure to match updated volume, elastance
        systems = map((a1,a2)->CVModule.elastancemodel!(a1,a2),systems,n);

        # output brachial pressure estimate
        Phat[numassims+1] = yhat[end]*Ps;
        println("Brachial pressure estimate: $(Phat[numassims+1])")
        numassims+=1;
    end
end
toc();

systems = pmap((a1,a2)->CVModule.updatevolumes!(a1,a2),systems,n);

# normalized RMSE ratio (optimal ensemble yields NRR ~ 1)
if assimflag == "yes"
    r2dot = sqrt.(1/systems[1].t[end]*r2dot);
    R1 = 1/systems[1].t[end]*r1dot;
    R2 = sum(r2dot);
    Ra = R1/R2;
    ERa = sqrt.((ensemblesize)/(2*(ensemblesize)));
    println("Normalized RMSE ratio: $(Ra/ERa)")
end

# if coupleflag == "yes"
#     ZMQ.close(sender)
#     ZMQ.close(ctx)
# end

if saveflag == "yes"
    fnames = ["p30_1_$i.mat" for i=1:ensemblesize];
    for i = 1:length(fnames)
        file = MAT.matopen(fnames[i],"w");
        write(file,"system",systems[i])
        close(file)
    end
end

if assimflag == "yes"
    return systems, n, βhat, τ1hat, Phat
else
    return systems, n
end

end
