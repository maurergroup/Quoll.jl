using Base

const AtomPairKeyedArray{T, N, AT, KT} = AbstractArray{<:KeyedArray{T, N, AT, KT}, 2}

abstract type AbstractCanonicalMetadata{A, S, B, P} <: AbstractOperatorMetadata{A, S, B, P} end
abstract type AbstractCanonicalOperator{O, T, D, M} <: AbstractOperator{O, T, D, M} end

abstract type AbstractBSparseMetadata{A, S, B, P} <: AbstractCanonicalMetadata{A, S, B, P} end
abstract type AbstractBSparseOperator{O, T, D, M, KD} <: AbstractCanonicalOperator{O, T, D, M} end

function Base.getindex(operator::AbstractBSparseOperator, key::NTuple{2, ChemicalSpecies})
    return getindex(operator.data, key)
end

function Base.setindex!(operator::AbstractBSparseOperator, val, key::NTuple{2, ChemicalSpecies})
    return getindex(operator.data, val, key)
end

function Base.getindex(operator::AbstractBSparseOperator, key::NTuple{2, Int})
    return getindex(operator.keydata, key)
end

function Base.setindex!(operator::AbstractBSparseOperator, val, key::NTuple{2, Int})
    return getindex(operator.keydata, val, key)
end

const CanonicalSHConversion = SHConversion(
    [[1], [1, 2, 3], [1, 2, 3, 4, 5], [1, 2, 3, 4, 5, 6, 7], [1, 2, 3, 4, 5, 6, 7, 8, 9]],
    [[1], [1, 1, 1], [1, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1, 1, 1, 1]],
)

SHConversion(::Type{<:AbstractCanonicalOperator}) = CanonicalSHConversion