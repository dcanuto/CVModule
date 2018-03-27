importall CVModule

function main()

rstflag = "no" # restart from prior solution, filename format must follow fnames below
hemoflag = "no" # 10% total blood vol. hemorrhage from left femoral artery
saveflag = "no" # save solution file to .mat struct
coupleflag = "no" # coupling with 3D organ model via ZMQ
assimflag = "yes" # patient data assimilation via EnKF to tune model state & params.
ptbflag = "yes" # generate ensemble via random perturbations, ONLY USE ONCE

ensemblesize = 3;
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

if assimflag == "yes"
    # allocators for ensemble augmented state, measurements
    X = [zeros(systems[1].solverparams.JL + length(systems[1].error.pdev)) for i in (1:length(systems)-1)];
    Y = [zeros(systems[1].solverparams.JL) for i in (1:length(systems)-1)];

    # allocators for state, measurement mean
    xhat = zeros(systems[1].solverparams.JL + length(systems[1].error.pdev));
    yhat = zeros(systems[1].solverparams.JL);
end

tic();
while systems[1].solverparams.numbeats < systems[1].solverparams.numbeatstotal
    soln = pmap((a1,a2)->CVModule.advancetime!(a1,a2;injury=hemoflag),systems,n);
    systems = [soln[i][1] for i in 1:length(soln)];
    n = [soln[i][2] for i in 1:length(soln)];
    if assimflag == "yes"
        # measurement
        y = systems[1].branches.P[61][n[1]+1,:];
        # ensemble forecast state and measurements
        X = [systems[i].branches.A[61][n[i]+1,:] for i in 2:length(soln)];
        for i = 2:length(soln)
            p1 = systems[i].branches.beta[61][end] + rand(Distributions.Normal(0,systems[1].error.pdev[1]));
            p2 = systems[i].branches.A0[61][end] + rand(Distributions.Normal(0,systems[1].error.pdev[2]));
            append!(X[i-1],[p1,p2])
        end
        println("Forecast X: $X")
        Y = [systems[i].branches.P[61][n[i]+1,:] for i in 2:length(soln)];
        println("Y: $Y")
        # ensemble mean state, measurement
        xhat = mean(X);
        yhat = mean(Y);
        println("̂x: $xhat")
        println("̂y: $yhat")
        # forecast state/meas. cross covariance, meas. covariance
        Pxy = zeros(length(xhat),length(yhat))
        Pyy = zeros(length(yhat),length(yhat))
        for i = 2:length(soln)
            Pxy += *((X[i-1] - xhat),(Y[i-1] - yhat)');
            Pyy += *((Y[i-1] - yhat),(Y[i-1] - yhat)');
        end
        Pxy /= (ensemblesize-1);
        Pyy /= (ensemblesize-1);
        # add meas. noise to meas. covariance (allows invertibility)
        Pyy += diagm(systems[1].error.odev[1]^2*ones(length(yhat)),0);
        # Kalman gain
        K = Pxy*inv(Pyy);
        # analysis step
        for i = 2:length(soln)
            X[i-1] += K*(y + rand(Distributions.Normal(0,systems[i-1].error.odev[1])) - Y[i-1]);
        end
        println("Analysis X: $X")
        # analysis back into ensemble members (TODO)
        println("Current time: $(systems[1].t[n[1]+1])")
    end
end
toc();

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
