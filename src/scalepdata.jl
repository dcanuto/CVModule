function scalepdata!(system::CVSystem,filename)
    file = MAT.matread(filename);
    pdata = file["pdata"];
    pdata = pdata[1,:];
    # println("Pressure data: $(pdata)")
    numassims = Int(cld(system.solverparams.numsteps,system.pdata.nsamp));
    pditp = Interpolations.interpolate(pdata, Interpolations.BSpline(Interpolations.Linear()), Interpolations.OnGrid());
    t = linspace(0,system.t[end],length(pdata));
    sitp = scale(pditp,t);
    tq = Float64[];
    for i = 1:numassims
        push!(tq,system.t[((i-1)*system.pdata.nsamp)+2])
    end
    # println("Final assimilation time: $(tq[end]) s")
    pq = [sitp[j] for j in tq];
    system.pdata.P = pq;
    return system
end
