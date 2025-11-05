@option struct BasisProjectionParams <: AbstractQuollParams
    projected_basis::Vector{BasisMetadata{Base.ImmutableDict{String, String}}}
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

    extras_pairs = (
        convert(String, match.captures[1]) => convert(String, match.captures[2])
        for match in eachmatch_extras
    )

    extras_dict_start = Base.ImmutableDict{String, String}()
    if isempty(extras_pairs)
        extras_dict = extras_dict_start
    else
        extras_dict = Base.ImmutableDict(extras_dict_start, extras_pairs...)
    end

    return (BasisMetadata(z, n, l, m, extras_dict) for m in ms)
end

function Configurations.from_dict(
    ::Type{BasisProjectionParams},
    ::OptionField{:projected_basis},
    ::Type{Vector{BasisMetadata{Base.ImmutableDict{String, String}}}},
    projected_basis_str,
)
    @argcheck isa(projected_basis_str, Vector{String}) "`projected_basis` field must be an array of strings"

    projected_basis = BasisMetadata{Base.ImmutableDict{String, String}}[]
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
    
    symbol = Symbol(normalize_comparison(s))
    @argcheck hasmethod(get_basis_projection, Tuple{Val{symbol}}) "Basis projection $s is unavailable or doesn't exist"

    return get_basis_projection(Val(symbol))
end
