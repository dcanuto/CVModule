type AssimErrors
    merr::Vector{Float64} # model error
    oerr::Vector{Float64} # obs. error
    perr::Vector{Float64} # param. error

    mdev::Vector{Float64} # σ of error distributions
    odev::Vector{Float64}
    pdev::Vector{Float64}

    mcv::Any # model covariance matrix
    ocv::Any # obs. covariance matrix

    mvar::Vector{Float64} # model variance (diag(mcv))
    ovar::Vector{Float64} # obs. variance

    a::Float64
    h::Float64

    lb::Vector{Float64} # lower/upper bounds for truncated normal param distrs.
    ub::Vector{Float64}

    pbar::Vector{Float64} # parameter distribution means

    function AssimErrors()
        this = new();
        this.merr = Vector{Float64}[];
        this.oerr = Vector{Float64}[];
        this.perr = Vector{Float64}[];
        this.mdev = Vector{Float64}[];
        this.odev = Vector{Float64}[];
        this.pdev = Vector{Float64}[];
        this.mvar = Vector{Float64}[];
        this.ovar = Vector{Float64}[];
        this.mcv = Array{Float64,2}[];
        this.ocv = Array{Float64,2}[];
        this.a = 1;
        this.h = 0;
        this.lb = Vector{Float64}[];
        this.ub = Vector{Float64}[];
        this.pbar = Vector{Float64}[];

        return this
    end
end
