abstract type AbstractCanonicalOperator <: AbstractOperator end
abstract type AbstractBSparseOperator <: AbstractCanonicalOperator end

const CanonicalSHConversion = SHConversion(
    [[1], [1, 2, 3], [1, 2, 3, 4, 5], [1, 2, 3, 4, 5, 6, 7], [1, 2, 3, 4, 5, 6, 7, 8, 9]],
    [[1], [1, 1, 1], [1, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1, 1, 1, 1]],
)

SHConversion(::Type{<:AbstractCanonicalOperator}) = CanonicalSHConversion