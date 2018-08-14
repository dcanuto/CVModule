function update0d!(system::CVSystem,n::Int64,terms::Vector{Int64},hemoflag="no")
    # update periphery
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

    # update heart chamber pressures
    CVModule.elastancemodel!(system,n+1);
end
