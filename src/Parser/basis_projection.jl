using Configurations
using AtomsBase
using ..Basis
using ..BasisProjection
using ..Utils

# TODO: Could make this parametric instead of forcing this into Dict{String, String}
# but for now parse_basismetadata() returns Dict{String, String} anyways
@option struct BasisProjectionParams <: AbstractQuollParams
    projected_basis::Vector{BasisMetadata{Dict{String, String}}}
    method::Type{<:AbstractBasisProjection}
end

# TODO: possibly implement as a BasisMetadata constructor?
function parse_basismetadata(basisfunc::AbstractString)
    re_znlm = r"(\w+)\s*\(\s*(\d+)\s*\,\s*(\d+)\s*\,\s*(\-?[\d+\*])\s*(\))"
    re_extras = r"(\w+)\s*\=\s*(\w+)"

    match_znlm = match(re_znlm, basisfunc)
    z = ChemicalSpecies(Symbol(match_znlm.captures[1]))
    n = parse(Int, match_znlm.captures[2])
    l = parse(Int, match_znlm.captures[3])

    m_str = match_znlm.captures[4]
    ms = m_str == "*" ? Tuple(-l:l) : (parse(Int, m_str),)

    extras_begin = match_znlm.offsets[5] + 1
    eachmatch_extras = eachmatch(re_extras, basisfunc[extras_begin:end])
    extras = Dict{String, String}(
        convert(String, match.captures[1]) => convert(String, match.captures[2])
        for match in eachmatch_extras
    )

    return (BasisMetadata(z, n, l, m, extras) for m in ms)
end

function Configurations.from_dict(
    ::Type{BasisProjectionParams},
    ::OptionField{:projected_basis},
    ::Type{Vector{BasisMetadata{Dict{String, String}}}},
    projected_basis_str,
)
    @argcheck isa(projected_basis_str, Vector{String}) "`projected_basis` field must be an array of strings"

    projected_basis = BasisMetadata{Dict{String, String}}[]
    for basisfunc in projected_basis_str
        push!(projected_basis, parse_basismetadata(basisfunc)...)
    end

    return projected_basis
end

function Configurations.from_dict(
    ::Type{BasisProjectionParams},
    ::OptionField{:method},
    ::Type{Type{<:AbstractBasisProjection}},
    s,
)
    @argcheck isa(s, String)
    
    symbol = Symbol(Utils.normalize_comparison(s))
    @argcheck hasmethod(basis_projection, Tuple{Val{symbol}}) "Basis projection $s is unavailable or doesn't exist"

    return basis_projection(Val(symbol))
end
