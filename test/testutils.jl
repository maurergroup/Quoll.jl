"""
    create_tempdirs(dirs)

Create a temporary tree of directories stored in a collection of
relative paths in `dirs`. Returns the root and the resulting paths
of the created directories.
"""
function create_tempdirs(dirs)
    root = mktempdir(; cleanup = true)

    paths = []
    for dir in dirs
        push!(paths, mkpath(joinpath(root, dir)))
    end

    return (root, paths)
end
