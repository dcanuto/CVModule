function sampletoparams!(system::CVSystem,params::Vector{Float64},restart::String,
    assim::String,sample::String)
    # arteries & periphery
    system.branches.beta[:] *= params[1];
    system.branches.A0[:] *= params[2];
    for i in 1:length(system.branches.ID)
        system.branches.c0[i] = [sqrt(0.5*system.branches.beta[i][end]/
            system.solverparams.rho)*system.branches.A0[i][end]^0.25];
    end
    Rdefault = [0.3,0.21,0.003,0.01]*mmHgToPa/cm3Tom3.*params[3:6];
    Cdefault = [0.01,1.64,1.81,13.24,73.88]*cm3Tom3/mmHgToPa.*params[7:11];
    Vdefault = [370.,370.,400.,500.,1400.]*cm3Tom3.*params[12:16];
    Ldefault = 5e-5*mmHgToPa/cm3Tom3*params[70];
    lowerflowfraction = 0.7*params[48];
    venousfractionofR = 0.1*params[49];
    venousfractionofL = 0.05*params[50];
    venousfractionofC = 0.95*params[51];
    venousfractionofV0 = 0.9*params[52];
    assignterminals!(system,Rdefault,Cdefault,Vdefault,Ldefault,lowerflowfraction,
        venousfractionofR,venousfractionofL,venousfractionofC,venousfractionofV0,
        Dict("a"=>0),restart,assim,sample);
    applybranchics!(system);
    applyperipheryics!(system);
    # liver
    system.liver.R[:] .*= params[17:20];
    system.liver.C[:] .*= params[21:24];
    system.liver.V0[:] .*= params[25:28];
    system.liver.L *= params[29];
    applyliverics!(system);
    # lungs
    system.lungs.Rp *= params[30];
    system.lungs.Ra .*= params[31:33];
    system.lungs.Rv .*= params[34:35];
    system.lungs.Ca .*= params[36:38];
    system.lungs.Cv .*= params[39:40];
    system.lungs.V0 .*= params[41:45];
    system.lungs.La *= params[46];
    system.lungs.Lv *= params[47];
    applylungics!(system);
    # heart
    system.heart.activation.m1 *= params[53];
    system.heart.activation.m2 *= params[54];
    system.heart.activation.tau1 *= params[55];
    system.heart.activation.tau2 *= params[56];
    t = linspace(0,system.heart.activation.th[1],100);
    g1 = (t/system.heart.activation.tau1).^system.heart.activation.m1;
    g2 = (t/system.heart.activation.tau2).^system.heart.activation.m2;
    system.heart.activation.k = [maximum((g1./(1+g1)).*(1./(1+g2)))^-1];
    system.heart.lv.Emax *= params[57];
    system.heart.rv.Emax *= params[58];
    system.heart.lv.V0 *= params[59];
    system.heart.rv.V0 *= params[60];
    system.heart.la.V0 *= params[61];
    system.heart.la.R *= params[62];
    system.heart.la.L *= params[63];
    system.heart.la.E *= params[64];
    system.heart.ra.V0 *= params[65];
    system.heart.ra.R *= params[66];
    system.heart.ra.L *= params[67];
    system.heart.ra.E *= params[68];
    system.heart.rv.L *= params[69];
    applyheartics!(system);
    # applycustomics!(system);
    updatevolumes!(system,0);

end
