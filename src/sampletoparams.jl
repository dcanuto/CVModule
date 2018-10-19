function sampletoparams(system::CVSystem,params::Vector{Float64},restart::String,
    assim::String,sample::String)
    # arteries & periphery
    system.branches.beta[:] *= params[1];
    system.brances.A0[:] *= params[2];
    Rdefault = [0.3,0.21,0.003,0.01]*mmHgToPa/cm3Tom3.*params[3:6];
    Cdefault = [0.01,1.64,1.81,13.24,73.88]*cm3Tom3/mmHgToPa.*params[7:10];
    Vdefault = [370.,370.,400.,500.,1400.]*cm3Tom3.*params[11:14];
    Ldefault = 5e-5*mmHgToPa/cm3Tom3*params[73];
    lowerflowfraction = 0.7*params[46];
    venousfractionofR = 0.1*params[47];
    venousfractionofL = 0.05*params[48];
    venousfractionofC = 0.95*params[49];
    venousfractionofV0 = 0.9*params[50];
    assignterminals!(system,Rdefault,Cdefault,Vdefault,Ldefault,lowerflowfraction,
        venousfractionofR,venousfractionofL,venousfractionofC,venousfractionofV0,
        Dict("a"=>0),restart,assim,sample);
    # liver
    system.liver.R[:] .*= params[15:18];
    system.liver.C[:] .*= params[19:22];
    system.liver.V0[:] .*= params[23:26];
    system.liver.L *= params[27];
    # lungs
    system.lungs.Rp *= params[28];
    system.lungs.Ra .*= params[29:31];
    system.lungs.Rv .*= params[32:33];
    system.lungs.Ca .*= params[34:36];
    system.lungs.Cv .*= params[37:38];
    system.lungs.V0 .*= params[39:43];
    system.lungs.La *= params[44];
    system.lungs.Lv *= params[45];
    # heart
    system.heart.activation.m1 *= params[51];
    system.heart.activation.m2 *= params[52];
    system.heart.activation.tau1 *= params[53];
    system.heart.activation.tau2 *= params[54];
    t = linspace(0,this.th[1],100);
    g1 = (t/this.tau1).^this.m1;
    g2 = (t/this.tau2).^this.m2;
    system.heart.activation.k = [maximum((g1./(1+g1)).*(1./(1+g2)))^-1];
    system.heart.lv.Emax *= params[55];
    system.heart.rv.Emax *= params[56];
    system.heart.lv.V0 *= params[57];
    system.heart.rv.V0 *= params[58];
    system.heart.la.V0 *= params[59];
    system.heart.la.R *= params[60];
    system.heart.la.L *= params[61];
    system.heart.la.E *= params[62];
    system.heart.ra.V0 *= params[63];
    system.heart.ra.R *= params[64];
    system.heart.ra.L *= params[65];
    system.heart.ra.E *= params[66];
    system.heart.rv.L *= params[67];
    system.heart.av.Kvo *= params[68];
    system.heart.av.Kvc *= params[69];
    system.heart.av.leff *= params[70];
    system.heart.av.Aann *= params[71];
    system.heart.av.Ks *= params[72];
    
end
