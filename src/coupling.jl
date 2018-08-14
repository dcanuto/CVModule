function coupling!(system::CVSystem,n::Int64,h::Float64,terms::Vector{Int64},hemoflag="no")
    if hemoflag == "no"
        # println("Distal coupling:")
        @time CVModule.coupledistal!(system,n,h,terms);
    elseif hemoflag == "yes"
        CVModule.coupledistal!(system,n,h,terms,hemoflag);
    end
    # println("Proximal coupling:")
    CVModule.coupleproximal!(system,n,h);
end
