using Configurations
using AtomsBase
using ..Basis
using ..BasisProjection
using ..Utils

@option struct BasisProjectionParams{E} <: AbstractQuollParams
    projected_basis::Vector{BasisMetadata{E}}
    method::AbstractBasisProjection
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
    ::Type{BasisProjectionParams{E}},
    ::OptionField{:projected_basis},
    ::Type{Vector{BasisMetadata{E}}},
    projected_basis_str,
) where E
    @argcheck isa(projected_basis_str, Vector{String}) "`projected_basis` field must be an array of strings"

    projected_basis = BasisMetadata[]
    for basisfunc in projected_basis_str
        push!(projected_basis, parse_basismetadata(basisfunc)...)
    end

    return Vector{typeof(first(projected_basis))}(projected_basis)
end

function Configurations.from_dict(
    ::Type{BasisProjectionParams{E}},
    ::OptionField{:method},
    ::Type{AbstractBasisProjection},
    method,
) where E
    @argcheck isa(method, String)

    itype = findfirst(
        x -> x == Utils.normalize_comparison(method),
        lowercase.(string.(nameof.(PROJ_REGISTRY)))
    )
    @argcheck itype !== nothing "Core projection method $method not found"

    return PROJ_REGISTRY[itype]()
end
