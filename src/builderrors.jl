type Errors
    a::Float64
    h::Float64

    odev::Vector{Float64}
    pdev::Vector{Float64}
    oscales::Vector{Float64}

    lb::Vector{Float64}
    ub::Vector{Float64}

    function Errors(nparams=1)
        this = new();
        δ = 0.99;
        this.a = (3*δ-1)/(2*δ);
        this.h = sqrt.(1-this.a^2);
        this.odev = [1e-12];
        this.pdev = 1e-12*ones(nparams);
        this.lb = -Inf*ones(nparams);
        this.ub = Inf*ones(nparams);
        this.oscales = ones(nparams);
        return this
    end
end
