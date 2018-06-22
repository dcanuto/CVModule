function regulateperiphery!(system::CVSystem,n::Int64,terms::Vector{Int64})
    for i = 1:length(terms)
        if system.branches.group[terms[i]] == "lower"
            system.branches.term[terms[i]].R[2] = system.cns.R2L[n+2];
            system.branches.term[terms[i]].R[3] = system.cns.R3L[n+2];
            system.branches.term[terms[i]].C[4] = system.cns.C4L[n+2];
            system.branches.term[terms[i]].C[5] = system.cns.C5L[n+2];
            system.branches.term[terms[i]].V0[4] = system.cns.V4L[n+2];
            system.branches.term[terms[i]].V0[5] = system.cns.V5L[n+2];
        else
            system.branches.term[terms[i]].R[2] = system.cns.R2U[n+2];
            system.branches.term[terms[i]].R[3] = system.cns.R3U[n+2];
            system.branches.term[terms[i]].C[4] = system.cns.C4U[n+2];
            system.branches.term[terms[i]].C[5] = system.cns.C5U[n+2];
            system.branches.term[terms[i]].V0[4] = system.cns.V4U[n+2];
            system.branches.term[terms[i]].V0[5] = system.cns.V5U[n+2];
        end
    end
end
