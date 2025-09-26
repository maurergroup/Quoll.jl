using Configurations
using ..Basis
using ..BasisProjection

@option struct BasisProjectionParams
    projected_basis::Vector{BasisMetadata}
    method::BasisProjection.AbstractBasisProjection
end

function Configurations.from_dict(
    ::Type{BasisProjectionParams},
    ::OptionField{:projected_basis},
    ::Vector{BasisMetadata},
    projected_basis,
)
    @argcheck isa(projected_basis, Vector{String}) "`projected_basis` field must be an array of strings"
    for basisfunc in projected_basis
        # TODO: possibly implement as a BasisMetadata constructor
        # "Si(2,    1,   -11)    extra   extra"
        re1 = r"\D{1,2}\s*\(\s*\d+\s*\,\s*\d+\s*\,\s*\-?[\d+\*]\s*\)"
        # Match the first part and deal with it, use split for the other one or match groups with \s+
    end
end

function Configurations.from_dict(
    ::Type{BasisProjectionParams},
    ::OptionField{:method},
    ::Type{BasisProjection.AbstractBasisProjection},
    method,
)
    itype = findfirst(x -> x == lowercase(method), lowercase.(string.(BasisProjection.IMPLEMENTED_PROJECTIONS)))
    @argcheck itype !== nothing "Core projection method $method not found"

    return BasisProjection.IMPLEMENTED_PROJECTIONS[itype]
end
