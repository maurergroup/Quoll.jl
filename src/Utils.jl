module Utils

"""
    find_leafdirs(rootdir::T) where T

Find leaf directories (directories that do not contain
directories themselves) going down from `rootdir` directory.
Absolute paths are returned.
"""
function find_leafdirs(rootdir::T) where {T}
    leafdirs = T[]
    for (root, dirs, _) in walkdir(rootdir)
        if isempty(dirs)
            push!(leafdirs, joinpath(rootdir, root))
        end
    end
    return leafdirs
end

end