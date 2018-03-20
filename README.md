# CVModule
Closed-loop model of the cardiovascular system with autonomic regulation. No fancy frills for interactive sessions yet (i.e., running a case requires manual changes to the `vascularmain.jl` script). To install:

    Pkg.clone("git://github.com/dcanuto/CVModule.git")

After installation, open the `vascularmain.jl` script and set options:

    rstflag = "yes" # (or "no")
    hemoflag = "yes"
    saveflag = "yes"
    
The above flags set whether or not to restart (in which case the `filename` variable must be set to the name of the restart file), to hemorrhage (set to perform a 10% hemorrhage and tourniquet application at the left femoral artery), and to save the solution to a restart file. For the last option, the restart file's name is set in the call to `matopen` near the end of the function, and can be edited as desired:

    file = MAT.matopen("savefile.mat", "w")

To run a case for `N` cardiac cycles, set the keyword argument `numbeatstotal` in the call to the solution initialization function within `vascularmain.jl`:

    system = CVModule.buildall(filename;numbeatstotal=N,restart=rstflag,injury=hemoflag);
    
With the above options set, build the `vascularmain.jl` script. To do so from the Julia REPL:

    include("path/to/vascularmain.jl")

Then run a case for `n` time steps from the REPL:

    system, n = main();
