function perturbics!(system::CVSystem)
    nptbs = 100; # number of perturbations
    fracmax = 0.1; # max. fraction of compartmental blood volume to switch
    termidxs = []; # indices of branches without children

    for i = 1:length(system.branches.ID)
        if isempty(system.branches.children[i])
            push!(termidxs,i);
        end
    end

    for i = 1:nptbs
        ptbidxs = shuffle(termidxs)[1:2];
        frac = rand()*fracmax;
        for j = 1:length(system.branches.term[ptbidxs[1]].V[1,:])
            system.branches.term[ptbidxs[2]].V[1,j] +=
                frac*system.branches.term[ptbidxs[1]].V[1,j];
            system.branches.term[ptbidxs[1]].V[1,j] *= (1-frac);
        end
    end

    return system
end
