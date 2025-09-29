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

"""
    normalize_comparison(s)

Normalize a string `s` by making the string lowercase and removing certain characters
(underscores and whitespace by default as defined in `replace_pairs`) for easier
comparison with reference values.
"""
function normalize_comparison(s::AbstractString; replace_pairs = ("_" => "", r"\s*" => ""))
    return lowercase(replace(s, replace_pairs...))
end

end