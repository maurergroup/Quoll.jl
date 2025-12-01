abstract type AbstractCanonicalMetadata{A, S, B, P} <: AbstractOperatorMetadata{A, S, B, P} end
abstract type AbstractCanonicalOperator{O, T, D, M} <: AbstractOperator{O, T, D, M} end

abstract type AbstractBSparseMetadata{A, S, B, P} <: AbstractCanonicalMetadata{A, S, B, P} end
abstract type AbstractDenseMetadata{A, S, B, P} <: AbstractCanonicalMetadata{A, S, B, P} end

const CanonicalSHConversion = SHConversion(
    [[1], [1, 2, 3], [1, 2, 3, 4, 5], [1, 2, 3, 4, 5, 6, 7], [1, 2, 3, 4, 5, 6, 7, 8, 9]],
    [[1], [1, 1, 1], [1, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1, 1, 1, 1]],
)

SHConversion(::Type{<:AbstractCanonicalOperator}) = CanonicalSHConversion

function Base.getindex(operator::AbstractCanonicalOperator, key::NTuple{2, ChemicalSpecies})
    return getindex(operator.data, key)
end

function Base.setindex!(operator::AbstractCanonicalOperator, val, key::NTuple{2, ChemicalSpecies})
    return getindex(operator.data, val, key)
end

function Base.getindex(operator::AbstractCanonicalOperator, key::NTuple{2, Int})
    return getindex(operator.keydata, key)
end

function Base.setindex!(operator::AbstractCanonicalOperator, val, key::NTuple{2, Int})
    return getindex(operator.keydata, val, key)
end

get_keydata(operator::AbstractCanonicalOperator) = operator.keydata

function construct_axis_metadata(z1::ChemicalSpecies, metadata::AbstractOperatorMetadata)
    
    basisset = get_basisset(metadata)
    spins = get_spins(metadata)

    if !isnothing(spins)
        axis_metadata = StructArray(
            orbital = basis_species(basisset, z1),
            spin = spins_species(spins, z1),
        )
    else
        axis_metadata = StructArray(
            orbital = basis_species(basisset, z1),
        )
    end
    
    return axis_metadata
end
