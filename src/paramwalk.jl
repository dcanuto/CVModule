function paramwalk!(system::CVSystem,θ::Vector{Float64})
    for i = 1:length(θ)
        θ[i] = rand(Distributions.TruncatedNormal(system.error.a*θ[i]+(1-system.error.a)*system.error.pbar[i],system.error.h*
            system.error.pdev[i],system.error.lb[i],system.error.ub[i]));
    end
    return θ
end
