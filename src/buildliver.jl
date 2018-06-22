type Liver # lumped portal vein, sinusoid, hepatic veins
    C::Vector{Float64}
    R::Vector{Float64}
    L::Float64
    V0::Vector{Float64}
    P::Any
    V::Any
    Q::Any

    function Liver()
        this = new()
        this.R = [0.004*mmHgToPa/cm3Tom3];
        append!(this.R,[0.005,0.005,0.004]*mmHgToPa/cm3Tom3);
        this.C = [3*cm3Tom3/mmHgToPa];
        append!(this.C,[10,15,45]*cm3Tom3/mmHgToPa);
        this.L = 5e-5*mmHgToPa/cm3Tom3;
        this.V0 = [50*cm3Tom3];
        append!(this.V0,[30,53,75]*cm3Tom3);
        this.P = Array{Float64,2}[];
        this.V = Array{Float64,2}[];
        this.Q = Array{Float64,2}[];
        return this
    end
end
