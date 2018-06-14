function applyendbcs!(system::CVSystem,n::Int64,terms::Vector{Int64},hemoflag="no")
    # update periphery
    if hemoflag == "no"
        CVModule.coupledistal!(system,n,terms);
    elseif hemoflag == "yes"
        CVModule.coupledistal!(system,n,terms,hemoflag);
    end
    for i = 1:length(terms)
        CVModule.updateterms!(system,n,terms[i]);
    end
    CVModule.updatevc!(system,n,terms);
    if hemoflag == "no"
        CVModule.updateliver!(system,n);
    elseif hemoflag == "yes"
        CVModule.updateliver!(system,n,hemoflag);
    end

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
