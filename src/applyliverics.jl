function applyliverics!(system::CVSystem,old=Dict("a"=>0),restart="no")
    if restart == "no"
        system.liver.P[1,1] = 4.99*mmHgToPa;
        system.liver.V[1,1] = (system.liver.P[1,1]*system.liver.C[1] +
            system.liver.V0[1]);
        system.liver.P[1,2] = 4.91*mmHgToPa;
        system.liver.V[1,2] = (system.liver.P[1,2]*system.liver.C[2] +
            system.liver.V0[2]);
        system.liver.P[1,3] = 4.79*mmHgToPa;
        system.liver.V[1,3] = (system.liver.P[1,3]*system.liver.C[3] +
            system.liver.V0[3]);
        system.liver.P[1,4] = 4.66*mmHgToPa;
        system.liver.V[1,4] = (system.liver.P[1,4]*system.liver.C[4] +
            system.liver.V0[4]);
    elseif restart == "yes"
        system.liver.P[1,1] = old["P"][end,1];
        system.liver.V[1,1] = old["V"][end,1];
        system.liver.Q[1,1] = old["Q"][end,1];
        system.liver.P[1,2] = old["P"][end,2];
        system.liver.V[1,2] = old["V"][end,2];
        system.liver.Q[1,2] = old["Q"][end,2];
        system.liver.P[1,3] = old["P"][end,3];
        system.liver.V[1,3] = old["V"][end,3];
        system.liver.Q[1,3] = old["Q"][end,3];
        system.liver.P[1,4] = old["P"][end,4];
        system.liver.V[1,4] = old["V"][end,4];
        system.liver.Q[1,4] = old["Q"][end,4];
    end
end
