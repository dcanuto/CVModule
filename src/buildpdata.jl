type PatientData
    dtsamp::Float64
    nsamp::Int64
    P::Vector{Float64}

    function PatientData()
        this = new();
        this.dtsamp = 1/125;
        nsamp = 0;
        P = Vector{Float64}[];
        return this
    end
end
