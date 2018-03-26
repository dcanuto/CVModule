type PatientData
    dtsamp::Float64
    nsamp::Int64
    P::Float64

    function PatientData()
        this = new();
        this.dtsamp = 0.1;
        nsamp = 0;
        P = 0;
        return this
    end
end
