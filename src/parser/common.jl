abstract type AbstractQuollParams end

function Configurations.from_dict(
    ::Type{OT},
    ::OptionField{:smearing_function},
    ::Type{Type{<:SmearingFunction}},
    s
) where OT<:AbstractQuollParams
    @argcheck isa(s, String)
    
    symbol = Symbol(normalize_comparison(s))
    @argcheck hasmethod(get_smearing, Tuple{Val{symbol}}) "Smearing function $s is unavailable or doesn't exist"

    return get_smearing(Val(symbol))
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
