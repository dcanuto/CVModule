function coupledistal!(system::CVSystem,n::Int64)
    # W1 at next time step
    CVModule.endinvariants!(system,n);

    for i = 1:length(system.branches.ID)
        if isempty(system.branches.children[i])
            CVModule.newtondist!(system,n,i);
        end
    end
end
