struct DeepHMetadata{A<:AbstractSystem, E} <: AbstractOperatorMetadata
    atoms::A
    sparsity::RealBlockSparsity
    basisset::BasisSetMetadata{E}
    spins::Union{SpinsMetadata, Nothing}
end

struct DeepHOperator{O<:AbstractOperatorKind, T<:AbstractFloat, A<:AbstractSystem, E} <: AbstractOperator
    kind::O
    data::Dictionary{NTuple{5, Int}, Array{T, 2}}
    metadata::DeepHMetadata{A, E}
end

# TODO: write
# - R_list.dat
# - element.dat
# - hamiltonians[_pred].h5
# - info.json
# - lat.dat
# - orbital_types.dat
# - rlat.dat
# - site_positions.dat

const DeepHSHConversion = SHConversion(
    [[1], [3, 1, 2], [3, 5, 1, 4, 2], [4, 5, 3, 6, 2, 7, 1]],
    [[1], [1, 1, 1], [1, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1, 1]],
)

SHConversion(::Type{DeepHOperator}) = DeepHSHConversion

get_writeformat(::Val{:deeph}) = DeepHOperator

get_avail_operatorkinds(::Type{DeepHOperator}) = [
    Hamiltonian(:ref, :none),
    Hamiltonian(:ref, :soc),
    Hamiltonian(:pred, :none),
    Hamiltonian(:pred, :soc),
    Overlap(:ref),
]

get_avail_filenames(::Hamiltonian, ::Val{:ref}, ::Val, ::Type{DeepHOperator}) = ["hamiltonians.h5"]
get_avail_filenames(::Hamiltonian, ::Val{:pred}, ::Val,  ::Type{DeepHOperator}) = ["hamiltonians_pred.h5"]

get_avail_filenames(::Overlap, ::Val{:ref}, ::Type{DeepHOperator}) = ["overlaps.h5"]

