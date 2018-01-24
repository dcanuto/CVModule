function loadtexttree(filename::String)
    branches = CSV.read(joinpath(pwd(),filename),null="");
    return branches
end
