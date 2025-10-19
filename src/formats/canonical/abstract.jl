using Base

abstract type AbstractCanonicalOperator <: AbstractOperator end
abstract type AbstractBSparseOperator <: AbstractCanonicalOperator end

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