# CardioModeling
Closed-loop model of the cardiovascular system with autonomic regulation. No fancy frills for interactive sessions yet (i.e., running a case requires manual changes to the `vascular_main.jl` function). To install:

    Pkg.clone("git://github.com/dcanuto/CardioModeling.git")

After installation, open the `vascular_main.jl` script and set options:

    rstflag = "yes" # (or "no")
    hemoflag = "yes"
    saveflag = "yes"
    
The above flags set whether or not to restart (in which case the `filename` variable must be set to the name of the restart file), to hemorrhage (set to perform a 10% hemorrhage and tourniquet application at the left femoral artery), and to save the solution. For the last option, the save file's name is set in the call to `matopen` near the end of the function:

    file = MAT.matopen("savefile.mat", "w")

To run a case for `N` cardiac cycles, set the keyword argument `numbeatstotal` in the call to the solution initialization function:

    system = CVModule.buildall(filename;numbeatstotal=N,restart=rstflag,injury=hemoflag);
    
With the above options set, build the `vascular_main` script. Then run a case for `n` time steps from the Julia REPL:

    system, n = main();
