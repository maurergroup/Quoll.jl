abstract type AbstractQuollParams end

function Configurations.from_dict(
    ::Type{OT},
    ::OptionField{:smearing_function},
    ::Type{Type{<:SmearingFunction}},
    s
) where OT<:AbstractQuollParams
    @argcheck isa(s, String)
    
    symbol = Symbol(Utils.normalize_comparison(s))
    @argcheck hasmethod(smearing, Tuple{Val{symbol}}) "Smearing function $s is unavailable or doesn't exist"

    return smearing(Val(symbol))
end

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

Normalize a string `s` by making the string lowercase and replacing certain characters
(underscores, hyphens and whitespace are removed by default as defined in `replace_pairs`)
for easier comparison with reference values.
"""
function normalize_comparison(s::AbstractString; replace_pairs = (rm => "" for rm in ["_", "-", r"\s*"]))
    return lowercase(replace(s, replace_pairs...))
end

