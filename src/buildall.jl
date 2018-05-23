type CVSystem # entire solution
    heart::CVModule.Heart
    branches::CVModule.ArterialBranches
    svc::CVModule.VenaCava
    ivc::CVModule.VenaCava
    lungs::CVModule.Lungs
    liver::CVModule.Liver
    hemo::CVModule.Hemorrhage
    cns::CVModule.CNS
    solverparams::CVModule.SolverParams
    t::Vector{Float64}
    arterialvolume::Float64
    peripheralvolume::Float64
    vcvolume::Float64
    heartvolume::Float64
    lungvolume::Float64
    livervolume::Float64
    initialvolume::Float64
    finalvolume::Float64

    function CVSystem(filename="test.csv",restart="no")
        this = new()
        if restart == "no"
            this.heart = CVModule.Heart();
            this.branches = CVModule.ArterialBranches(filename);
        elseif restart == "yes"
            vars = MAT.matread(filename);
            sys = vars["system"];
            heart = sys["heart"];
            branches = sys["branches"];
            this.heart = CVModule.Heart(heart,restart)
            this.branches = CVModule.ArterialBranches(filename,branches,restart)
        else
            error("Keyword restart must either be yes or no. Aborting.")
        end
        this.svc = CVModule.VenaCava();
        this.ivc = CVModule.VenaCava();
        this.lungs = CVModule.Lungs();
        this.liver = CVModule.Liver();
        this.cns = CVModule.CNS();
        this.hemo = CVModule.Hemorrhage();
        this.solverparams = CVModule.SolverParams();
        this.t = Vector{Float64}[];
        return this
    end
end

# build solution struct
function buildall(filename="test.csv";numbeatstotal=1,restart="no",injury="no")
    if restart == "no"
        system = CVModule.CVSystem(filename);
        system.solverparams.numbeatstotal = numbeatstotal;
        CVModule.calcbranchprops!(system);
        CVModule.discretizebranches!(system);
        CVModule.assignterminals!(system);
        CVModule.discretizeperiphery!(system);
        CVModule.discretizeheart!(system);
        CVModule.discretizelungs!(system);
        CVModule.discretizeliver!(system);
        CVModule.discretizecns!(system);
        CVModule.applybranchics!(system);
        CVModule.applyperipheryics!(system);
        CVModule.applyheartics!(system);
        CVModule.applylungics!(system);
        CVModule.applyliverics!(system);
        CVModule.applycnsics!(system);
        CVModule.applycustomics!(system);
        CVModule.applyhemoics!(system);
    elseif restart == "yes"
        vars = MAT.matread(filename);
        sys = vars["system"];
        branches = sys["branches"];
        term = branches["term"];
        heart = sys["heart"];
        lungs = sys["lungs"];
        liver = sys["liver"];
        cns = sys["cns"];
        system = CVModule.CVSystem(filename,restart);
        system.solverparams.numbeatstotal = numbeatstotal;
        CVModule.calcbranchprops!(system,branches,restart);
        CVModule.discretizebranches!(system,sys,restart);
        CVModule.assignterminals!(system,term,restart);
        CVModule.discretizeperiphery!(system);
        CVModule.discretizeheart!(system);
        CVModule.discretizelungs!(system);
        CVModule.discretizeliver!(system);
        CVModule.discretizecns!(system);
        CVModule.applybranchics!(system,sys,restart);
        CVModule.applyperipheryics!(system,sys,restart);
        CVModule.applyheartics!(system,heart,restart);
        CVModule.applylungics!(system,lungs,restart);
        CVModule.applyliverics!(system,liver,restart);
        CVModule.applycnsics!(system,cns,restart);
        CVModule.applyhemoics!(system,sys);
    end
    CVModule.updatevolumes!(system,0);
    return system
end
