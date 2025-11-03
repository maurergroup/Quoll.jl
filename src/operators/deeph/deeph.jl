const DeepHDict{T, N} = AbstractAnyDict{NTuple{5, Int}, <:AbstractArray{T, N}}

struct DeepHMetadata{A<:AbstractSystem, S<:RealBlockSparsity, B<:BasisSetMetadata, P<:Union{SpinsMetadata, Nothing}} <: AbstractOperatorMetadata{A, S, B, P}
    atoms::A
    sparsity::S
    basisset::B
    spins::P
end

struct DeepHOperator{O<:OperatorKind, T<:Number, D<:DeepHDict{T}, M<:DeepHMetadata} <: AbstractOperator{O, T, D, M}
    kind::O
    data::D
    metadata::M
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
    Hamiltonian(source = :ref, spin = :none),
    Hamiltonian(source = :ref, spin = :soc),
    Hamiltonian(source = :pred, spin = :none),
    Hamiltonian(source = :pred, spin = :soc),
    Overlap(source = :ref),
]

get_avail_filenames(
    ::Hamiltonian,
    ::Pair{Val{:source}, Val{:ref}},
    ::Pair{Val{:spin}, Val{SPIN}},
    ::Type{DeepHOperator}
) where SPIN = ["hamiltonians.h5"]

get_avail_filenames(
    ::Hamiltonian,
    ::Pair{Val{:source}, Val{:pred}},
    ::Pair{Val{:spin}, Val{SPIN}},
    ::Type{DeepHOperator}
) where SPIN = ["hamiltonians_pred.h5"]

get_avail_filenames(
    ::Overlap, 
    ::Pair{Val{:source}, Val{:ref}},
    ::Type{DeepHOperator}
) = ["overlaps.h5"]
