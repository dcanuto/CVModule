function applyendbcs!(system::CVSystem,n::Int64,hemoflag="no")
    # update periphery
    if hemoflag == "no"
        CVModule.coupledistal!(system,n);
    elseif hemoflag == "yes"
        CVModule.coupledistal!(system,n,hemoflag);
    end  
    CVModule.updateterms!(system,n);
    CVModule.updatevc!(system,n);

    # update right heart
    CVModule.updaterh!(system,n);

    # update pulmonary circulation
    CVModule.updatelungs!(system,n);

    # update left heart
    CVModule.updatela!(system,n);
    CVModule.coupleproximal!(system,n);

    # update heart chamber pressures
    CVModule.elastancemodel!(system,n+1);
end
