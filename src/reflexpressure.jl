function reflexpressure!(system::CVSystem,n::Int64)
    # ending time step index for averaging
    system.cns.avgindexend = n+2;

    # mean of internal carotid/aortic arch means
    avg1 = 1/(system.t[system.cns.avgindexend]-
        system.t[system.cns.avgindexstart])*NumericalIntegration.integrate(
        system.t[system.cns.avgindexstart:system.cns.avgindexend],
        system.branches.P[2][system.cns.avgindexstart:system.cns.avgindexend,1]
        );
    avg2 = 1/(system.t[system.cns.avgindexend]-
        system.t[system.cns.avgindexstart])*NumericalIntegration.integrate(
        system.t[system.cns.avgindexstart:system.cns.avgindexend],
        system.branches.P[53][system.cns.avgindexstart:system.cns.avgindexend,1]
        );
    avg3 = 1/(system.t[system.cns.avgindexend]-
        system.t[system.cns.avgindexstart])*NumericalIntegration.integrate(
        system.t[system.cns.avgindexstart:system.cns.avgindexend],
        system.branches.P[89][system.cns.avgindexstart:system.cns.avgindexend,1]
        );
    push!(system.cns.Paverage,round(1/3*(avg1 + avg2 + avg3),1)*mmHgToPa);
    println("Average baroreflex pressure over cardiac cycle: $(system.cns.Paverage[end]/mmHgToPa) mmHg")

    # beginning time step index for averaging
    system.cns.avgindexstart = system.cns.avgindexend;
end
