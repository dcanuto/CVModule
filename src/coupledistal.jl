function coupledistal!(system::CVSystem,n::Int64,hemoflag="no")
    # W1 at next time step
    CVModule.endinvariants!(system,n);
    if hemoflag == "no"
        for i = 1:length(system.branches.ID)
            if isempty(system.branches.children[i])
                CVModule.newtondist!(system,n,i);
            end
        end
    elseif hemoflag == "yes"
        for i = 1:length(system.branches.ID)
            if isempty(system.branches.children[i]) && system.hemo.injured[i] == false
                CVModule.newtondist!(system,n,i);
            elseif isempty(system.branches.children[i]) && system.hemo.injured[i] == true
                CVModule.model1dhemo!(system,n,i);
            end
        end
    end
end
