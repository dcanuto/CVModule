function updateterms!(system::CVSystem,n::Int64)
    for i = 1:length(system.branches.ID)
        if isempty(system.branches.children[i])
            # arteriolar compartment V, P
            system.branches.term[i].V[n+2,2] = (
                system.branches.term[i].V[n+1,2]+
                system.solverparams.h*(system.branches.term[i].Q[n+1,1]-
                system.branches.term[i].Q[n+1,2]));
            system.branches.term[i].P[n+2,2] = (
                1/system.branches.term[i].C[2]*(
                system.branches.term[i].V[n+2,2]-
                system.branches.term[i].V0[2]));
            # capillary compartment V, P
            system.branches.term[i].V[n+2,3] = (
                system.branches.term[i].V[n+1,3]+
                system.solverparams.h*(system.branches.term[i].Q[n+1,2]-
                system.branches.term[i].Q[n+1,3]));
            system.branches.term[i].P[n+2,3] = (
                1/system.branches.term[i].C[3]*(
                system.branches.term[i].V[n+2,3]-
                system.branches.term[i].V0[3]));
            # venular compartment V, P
            system.branches.term[i].V[n+2,4] = (
                system.branches.term[i].V[n+1,4]+
                system.solverparams.h*(system.branches.term[i].Q[n+1,3]-
                system.branches.term[i].Q[n+1,4]));
            system.branches.term[i].P[n+2,4] = (
                1/system.branches.term[i].C[4]*(
                system.branches.term[i].V[n+2,4]-
                system.branches.term[i].V0[4]));
            # venous compartment
            if i == 12 # inferior mesenteric flows into splenic
                system.branches.term[i].V[n+2,5] = (
                    system.branches.term[i].V[n+1,5] + system.solverparams.h*
                    (system.branches.term[21].Q[n+1,5] + system.branches.term[i].Q[n+1,4]-
                    system.branches.term[i].Q[n+1,5]));
            else
                system.branches.term[i].V[n+2,5] = (
                    system.branches.term[i].V[n+1,5]+
                    system.solverparams.h*(system.branches.term[i].Q[n+1,4]-
                    system.branches.term[i].Q[n+1,5]));
            end
            system.branches.term[i].P[n+2,5] = (
                1/system.branches.term[i].C[5]*(
                system.branches.term[i].V[n+2,5]-
                system.branches.term[i].V0[5]));
            if system.branches.group[i] == "lower"
                if i == 12 || i == 13 || i == 15 # divert flow to hepatic portal vein
                    system.branches.term[i].Q[n+2,5] = (
                        (1 - system.branches.term[i].R[5]/
                        system.branches.term[i].L[5]*system.solverparams.h)*
                        system.branches.term[i].Q[n+1,5]+
                        system.solverparams.h/system.branches.term[i].L[5]*(
                        system.branches.term[i].P[n+1,5]-
                        system.liver.P[n+1,1]));
                elseif i == 21 # inferior mesenteric joins splenic
                    system.branches.term[i].Q[n+2,5] = (
                        (1 - system.branches.term[i].R[5]/
                        system.branches.term[i].L[5]*system.solverparams.h)*
                        system.branches.term[i].Q[n+1,5]+
                        system.solverparams.h/system.branches.term[i].L[5]*(
                        system.branches.term[i].P[n+1,5]-
                        system.branches.term[12].P[n+1,5]));
                else
                    system.branches.term[i].Q[n+2,5] = (
                        (1 - system.branches.term[i].R[5]/
                        system.branches.term[i].L[5]*system.solverparams.h)*
                        system.branches.term[i].Q[n+1,5]+
                        system.solverparams.h/system.branches.term[i].L[5]*(
                        system.branches.term[i].P[n+1,5]-
                        system.ivc.P[n+1]));
                end
            else
                system.branches.term[i].Q[n+2,5] = (
                    (1 - system.branches.term[i].R[5]/
                    system.branches.term[i].L[5]*system.solverparams.h)*
                    system.branches.term[i].Q[n+1,5]+
                    system.solverparams.h/system.branches.term[i].L[5]*(
                    system.branches.term[i].P[n+1,5]-
                    system.svc.P[n+1]));
            end
            # Q for other compartments
            system.branches.term[i].Q[n+2,1] = (
                (system.branches.term[i].P[n+2,1]-
                system.branches.term[i].P[n+2,2])/
                system.branches.term[i].R[1]);
            system.branches.term[i].Q[n+2,2] = (
                (system.branches.term[i].P[n+2,2]-
                system.branches.term[i].P[n+2,3])/
                system.branches.term[i].R[2]);
            system.branches.term[i].Q[n+2,3] = (
                (system.branches.term[i].P[n+2,3]-
                system.branches.term[i].P[n+2,4])/
                system.branches.term[i].R[3]);
            system.branches.term[i].Q[n+2,4] = (
                (system.branches.term[i].P[n+2,4]-
                system.branches.term[i].P[n+2,5])/
                system.branches.term[i].R[4]);
        end
    end
end
