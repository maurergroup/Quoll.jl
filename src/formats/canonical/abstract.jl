abstract type AbstractCanonicalOperator <: AbstractOperator end
abstract type AbstractBSparseOperator <: AbstractCanonicalOperator end

struct BSparseOperatorMetadata <: AbstractOperatorMetadata
    z1z2_ij2interval::Dictionary{NTuple{2, ChemicalSpecies}, Dictionary{NTuple{2, Int64}, UnitRange{Int64}}}
    atom2species::Vector{ChemicalSpecies}
    # TODO: I think I can keep those three as separate structs and only keep
    # attributes which are required for matrix operations in this struct
    # Sparsity?
    # Atoms?
    # BasisSetMetadata?
    # TODO: possibly construct this by only taking BasisSetMetadata, Atoms, Sparsity,
    # but create an outer constructor which could compute z1z2_ij2interval from those internally
    # TODO: also possibly create a constructor which could create a matrix with zero values
end

struct RealBSparseOperator{O<:AbstractOperatorKind, T<:AbstractFloat, AT, KT} <: AbstractBSparseOperator
    kind::O
    data::Dictionary{NTuple{2, ChemicalSpecies}, Array{T, 3}}
    keydata::Dictionary{NTuple{2, Int}, KeyedArray{T, 3, AT, KT}}
    metadata::BSparseOperatorMetadata
end