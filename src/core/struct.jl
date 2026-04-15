### TRAITS ###

"""Base type for all Holy Traits used in dispatch throughout Quoll."""
abstract type AbstractTrait end

### METADATA TRAITS ###

"""Trait distinguishing real-space and reciprocal-space metadata."""
abstract type SpaceTrait <: AbstractTrait end

"""Singleton indicating real-space metadata (no k-point)."""
struct RealSpace <: SpaceTrait end

"""Singleton indicating reciprocal-space metadata (has a k-point)."""
struct RecipSpace <: SpaceTrait end

"""Trait distinguishing spin-polarised and non-spin-polarised metadata."""
abstract type SpinTrait <: AbstractTrait end

"""Singleton indicating no spin polarisation."""
struct NoSpin <: SpinTrait end

"""Singleton indicating spin-polarised metadata (carries a `SpinsMetadata`)."""
struct HasSpin <: SpinTrait end

### BASIC METADATA CONTAINER ###

"""
    BasicMetadataContainer{O,X,S,B,Y,A}

Holds the core metadata fields shared by all metadata types: operator kind, source format,
sparsity pattern, basis set, spherical harmonics convention, and atomic structure.

Every concrete `AbstractMetadata` subtype wraps a `BasicMetadataContainer` in its `common`
field, accessed via [`op_basic_metadata`](@ref).
"""
struct BasicMetadataContainer{
    O<:OperatorKind,
    X<:AbstractSource,
    S<:AbstractSparsity,
    B<:BasisSetMetadata,
    Y<:SHConvention,
    A<:AbstractSystem,
}
    kind::O
    source::X
    sparsity::S
    basisset::B
    shconv::Y
    atoms::A
end

op_kind(common::BasicMetadataContainer) = common.kind
op_source(common::BasicMetadataContainer) = common.source
op_sparsity(common::BasicMetadataContainer) = common.sparsity
op_basisset(common::BasicMetadataContainer) = common.basisset
op_shconv(common::BasicMetadataContainer) = common.shconv
op_atoms(common::BasicMetadataContainer) = common.atoms

op_hermicity(common::BasicMetadataContainer) = op_hermicity(op_sparsity(common))

### METADATA TYPES ###

"""
    AbstractMetadata{O,X,S,B,Y,A}

Abstract supertype for all operator metadata. Type parameters encode the operator kind (`O`),
source format (`X`), sparsity pattern (`S`), basis set (`B`), SH convention (`Y`), and
atomic structure (`A`).

Concrete subtypes, such as [`RealMetadata`](@ref), [`RecipMetadata`](@ref),
[`SpinRealMetadata`](@ref), [`SpinRecipMetadata`](@ref), differ in which extra fields
they carry (k-point, spin information, or both). Traits, e.g. `SpaceTrait` and `SpinTrait`,
encode these differences for dispatch.
"""
abstract type AbstractMetadata{O,X,S,B,Y,A} end

"""
    RealMetadata{O,X,S,B,Y,A}

Real-space, non-spin-polarised metadata. Wraps only a `BasicMetadataContainer`.
"""
struct RealMetadata{
    O<:OperatorKind,
    X<:AbstractSource,
    S<:AbstractSparsity,
    B<:BasisSetMetadata,
    Y<:SHConvention,
    A<:AbstractSystem,
} <: AbstractMetadata{O,X,S,B,Y,A}
    common::BasicMetadataContainer{O,X,S,B,Y,A}
end
trait(::Type{SpaceTrait}, ::Type{<:RealMetadata}) = RealSpace()
trait(::Type{SpinTrait}, ::Type{<:RealMetadata}) = NoSpin()
extrafield_traittypes(::Type{<:RealMetadata}) = []

"""
    RecipMetadata{O,X,S,B,Y,A,K}

Reciprocal-space, non-spin-polarised metadata. Extends `BasicMetadataContainer` with a
k-point vector.
"""
struct RecipMetadata{
    O<:OperatorKind,
    X<:AbstractSource,
    S<:AbstractSparsity,
    B<:BasisSetMetadata,
    Y<:SHConvention,
    A<:AbstractSystem,
    K<:SVector{3},
} <: AbstractMetadata{O,X,S,B,Y,A}
    common::BasicMetadataContainer{O,X,S,B,Y,A}
    kpoint::K
end
trait(::Type{SpaceTrait}, ::Type{<:RecipMetadata}) = RecipSpace()
trait(::Type{SpinTrait}, ::Type{<:RecipMetadata}) = NoSpin()
extrafield_traittypes(::Type{<:RecipMetadata}) = [SpaceTrait]

"""
    SpinRealMetadata{O,X,S,B,Y,A,P}

Real-space, spin-polarised metadata. Extends `BasicMetadataContainer` with spin information.
"""
struct SpinRealMetadata{
    O<:OperatorKind,
    X<:AbstractSource,
    S<:AbstractSparsity,
    B<:BasisSetMetadata,
    Y<:SHConvention,
    A<:AbstractSystem,
    P<:SpinsMetadata,
} <: AbstractMetadata{O,X,S,B,Y,A}
    common::BasicMetadataContainer{O,X,S,B,Y,A}
    spins::P
end
trait(::Type{SpaceTrait}, ::Type{<:SpinRealMetadata}) = RealSpace()
trait(::Type{SpinTrait}, ::Type{<:SpinRealMetadata}) = HasSpin()
extrafield_traittypes(::Type{<:SpinRealMetadata}) = [SpinTrait]

"""
    SpinRecipMetadata{O,X,S,B,Y,A,P,K}

Reciprocal-space, spin-polarised metadata. Extends `BasicMetadataContainer` with both
spin information and a k-point vector.
"""
struct SpinRecipMetadata{
    O<:OperatorKind,
    X<:AbstractSource,
    S<:AbstractSparsity,
    B<:BasisSetMetadata,
    Y<:SHConvention,
    A<:AbstractSystem,
    P<:SpinsMetadata,
    K<:SVector{3},
} <: AbstractMetadata{O,X,S,B,Y,A}
    common::BasicMetadataContainer{O,X,S,B,Y,A}
    spins::P
    kpoint::K
end
trait(::Type{SpaceTrait}, ::Type{<:SpinRecipMetadata}) = RecipSpace()
trait(::Type{SpinTrait}, ::Type{<:SpinRecipMetadata}) = HasSpin()
extrafield_traittypes(::Type{<:SpinRecipMetadata}) = [SpaceTrait, SpinTrait]

# Methods for optional larger unions
extrafield_traittypes(::Type{<:Union{<:RealMetadata,<:SpinRealMetadata}}) = [SpinTrait]
extrafield_traittypes(::Type{<:Union{<:RecipMetadata,<:SpinRecipMetadata}}) =
    [SpaceTrait, SpinTrait]

function Base.show(io::IO, metadata::AbstractMetadata)
    print(io, nameof(typeof(metadata)))
    print(io, "($(op_kind(metadata)), ")
    print(io, "$(nameof(typeof(op_source(metadata)))), ")
    return print(io, "$(nameof(typeof(op_sparsity(metadata)))))")
end

### METADATA ACCESSORS ###

op_basic_metadata(metadata::AbstractMetadata) = metadata.common

op_kpoint(metadata::AbstractMetadata) = metadata.kpoint
op_spins(metadata::AbstractMetadata) = metadata.spins

op_kind(metadata::AbstractMetadata) = op_kind(op_basic_metadata(metadata))
op_source(metadata::AbstractMetadata) = op_source(op_basic_metadata(metadata))
op_sparsity(metadata::AbstractMetadata) = op_sparsity(op_basic_metadata(metadata))
op_basisset(metadata::AbstractMetadata) = op_basisset(op_basic_metadata(metadata))
op_shconv(metadata::AbstractMetadata) = op_shconv(op_basic_metadata(metadata))
op_atoms(metadata::AbstractMetadata) = op_atoms(op_basic_metadata(metadata))

op_hermicity(metadata::AbstractMetadata) = op_hermicity(op_basic_metadata(metadata))

### DATA ###

"""
    DataContainer{T,N,B,X,S}

Type-tagged wrapper around operator data arrays. The type parameters encode the element type
(`T`), dimensionality (`N`), concrete body type (`B`), source format (`X`), and sparsity
pattern (`S`). This enables dispatch on data format without inspecting the body directly.

Use [`unwrap_data`](@ref) to access the underlying array/dictionary and
[`wrap_data`](@ref) to construct a `DataContainer` from a raw body.
"""
struct DataContainer{T,N,B,X<:AbstractSource,S<:AbstractSparsity}
    body::B
end

"""
    unwrap_data(data::DataContainer) -> body

Return the underlying array or dictionary stored in the `DataContainer`.
"""
unwrap_data(data::DataContainer) = data.body

### COMMON DATA TYPES ###

const CSCRealData{T} = DataContainer{T,1,<:AbstractVector{T},<:Any,<:CSCRealSparsity}
const DenseRecipData{T} = DataContainer{T,2,<:AbstractMatrix{T},<:Any,<:DenseRecipSparsity}

"""
    wrap_data(::Type{M}, body) where {M<:AbstractMetadata}
    wrap_data(::Type{X}, ::Type{S}, body)

Wrap a raw array or dictionary into a [`DataContainer`](@ref), tagging it with the source
and sparsity types inferred from metadata type `M` (or provided directly as `X` and `S`).
"""
function wrap_data(::Type{M}, body) where {M<:AbstractMetadata}
    return wrap_data(op_source_type(M), op_sparsity_type(M), body)
end

function wrap_data(
    ::Type{X}, ::Type{S}, body::B
) where {
    X<:AbstractSource,
    S<:AbstractSparsity,
    B<:AbstractAnyDict{K,V},
} where {K,V<:AbstractArray{T,N}} where {T<:Number,N}
    return DataContainer{T,N,B,X,S}(body)
end

function wrap_data(
    ::Type{X}, ::Type{S}, body::B
) where {
    X<:AbstractSource,
    S<:AbstractSparsity,
    B<:AbstractArray{T,N},
} where {T<:Number,N}
    return DataContainer{T,N,B,X,S}(body)
end

### OPERATOR TRAITS ###

"""Trait distinguishing operators that carry keydata from those that do not."""
abstract type KeyedTrait <: AbstractTrait end

"""Singleton indicating an operator carries keyed data (e.g. per-atom-pair blocks)."""
struct HasKeydata <: KeyedTrait end

"""Singleton indicating an operator has no keyed data."""
struct NoKeydata <: KeyedTrait end

### OPERATOR ###

"""
    AbstractOperator{O,M,D}

Abstract supertype for all operators. Type parameters encode the operator kind (`O`),
metadata type (`M`), and data container type (`D`).
"""
abstract type AbstractOperator{O,M,D} end

"""
    Operator{O,M,D}

An operator that holds metadata and a single data container.

`trait(KeyedTrait, Operator) → NoKeydata()`.

See also [`KeyedOperator`](@ref), [`build_operator`](@ref).
"""
struct Operator{
    O<:OperatorKind,
    M<:AbstractMetadata{O},
    D<:DataContainer,
} <: AbstractOperator{O,M,D}
    metadata::M
    data::D

    function Operator(
        metadata::M, data::D
    ) where {M<:AbstractMetadata{O},D<:DataContainer} where {O}
        return new{O,M,D}(metadata, data)
    end
end
trait(::Type{KeyedTrait}, ::Type{<:Operator}) = NoKeydata()
extrafield_traittypes(::Type{<:Operator}) = []

"""
    KeyedOperator{O,M,D,KD}

An operator with keyed data. Holds metadata, a data container, and an additional keydata
container (e.g. per-atom-pair block matrices indexed by `(i, j)` pairs).

`trait(KeyedTrait, KeyedOperator) → HasKeydata()`.

See also [`Operator`](@ref), [`build_operator`](@ref).
"""
struct KeyedOperator{
    O<:OperatorKind,
    M<:AbstractMetadata{O},
    D<:DataContainer,
    KD<:DataContainer,
} <: AbstractOperator{O,M,D}
    metadata::M
    data::D
    keydata::KD

    function KeyedOperator(
        metadata::M, data::D, keydata::KD
    ) where {M<:AbstractMetadata{O},D<:DataContainer,KD<:DataContainer} where {O}
        return new{O,M,D,KD}(metadata, data, keydata)
    end
end
trait(::Type{KeyedTrait}, ::Type{<:KeyedOperator}) = HasKeydata()
extrafield_traittypes(::Type{<:KeyedOperator}) = [KeyedTrait]

function Base.show(io::IO, operator::AbstractOperator)
    print(io, nameof(typeof(operator)))
    print(io, "(")
    show(io, op_metadata(operator))
    return print(io, ")")
end

"""
    synchronise_data!(operator, comm::MPI.Comm)

All-reduce the operator's data across MPI communicator `comm`, summing contributions from
different ranks (e.g. after distributing k-points across processes).
"""
function synchronise_data!(operator::AbstractOperator, comm::MPI.Comm)
    data = op_data(operator)
    return synchronise_data!(data, comm)
end

### OPERATOR ACCESSORS ###

op_metadata_type(::AbstractOperator{O,M}) where {O,M} = M
op_data_type(::AbstractOperator{O,M,D}) where {O,M,D} = D
op_float_type(::AbstractOperator{O,M,D}) where {O,M,D<:DataContainer{T}} where {T} = T

op_data(operator::AbstractOperator) = operator.data
op_metadata(operator::AbstractOperator) = operator.metadata

op_keydata(operator::AbstractOperator) = operator.keydata

op_kind(operator::AbstractOperator) = op_kind(op_metadata(operator))
op_source(operator::AbstractOperator) = op_source(op_metadata(operator))
op_sparsity(operator::AbstractOperator) = op_sparsity(op_metadata(operator))
op_basisset(operator::AbstractOperator) = op_basisset(op_metadata(operator))
op_shconv(operator::AbstractOperator) = op_shconv(op_metadata(operator))
op_atoms(operator::AbstractOperator) = op_atoms(op_metadata(operator))

op_hermicity(operator::AbstractOperator) = op_hermicity(op_metadata(operator))

op_kpoint(operator::AbstractOperator) = op_kpoint(op_metadata(operator))
op_spins(operator::AbstractOperator) = op_spins(op_metadata(operator))