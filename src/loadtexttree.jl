function loadtexttree(filename::String)
    branches = CSV.read(joinpath(abspath(dirname(@__FILE__)),filename),null="");
    return branches
end
