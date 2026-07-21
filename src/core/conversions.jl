### OPERATOR ###

"""
    convert_operator(::Type{OPₒᵤₜ}, ::Type{Mₒᵤₜ}, in_operator; kwargs...) -> OPₒᵤₜ

Convert an operator to a different format in a single call. This is the main entry point
for format conversion, composing three stages:

1. [`convert_metadata`](@ref) — convert source, SH convention, sparsity, basis set.
2. [`build_operator`](@ref) — allocate the output operator with zero-initialised data.
3. [`convert_data!`](@ref) — transfer and transform the actual matrix data.

# Arguments
- `OPₒᵤₜ`: output operator type (`<:KeyedOperator`).
- `Mₒᵤₜ`: output metadata type (e.g. `CanonicalBlockRealMetadata`, `DeepHBlockRealMetadata`).
  May be a union type; the concrete subtype is chosen based on available extra fields.
- `in_operator`: the input operator to convert.

# Keyword arguments
- `radii=nothing`: if given, build output sparsity from a neighbour list with these radii
  instead of converting the input sparsity.
- `hermitian=nothing`: override hermicity of the output (defaults to input's hermicity).
- `float=nothing`: override the floating-point element type (defaults to input's).
- `out_shconv=nothing`: explicit output SH convention (defaults to the output source's default).
- `source_kwargs=NamedTuple()`: extra keyword arguments for constructing the output source struct.
- `operator_extra_kwargs=NamedTuple()`: passed to `build_operator_extra`.
- `metadata_extra_kwargs=NamedTuple()`: passed to `convert_metadata_extra` (e.g. `kpoint`).
"""
function convert_operator(
    ::Type{OPₒᵤₜ}, ::Type{Mₒᵤₜ}, in_operator::OPᵢₙ;
    radii=nothing, hermitian=nothing, float=nothing, out_shconv=nothing,
    source_kwargs=NamedTuple(), operator_extra_kwargs=NamedTuple(),
    metadata_extra_kwargs=NamedTuple(),
) where {OPₒᵤₜ<:AbstractOperator,Mₒᵤₜ<:AbstractMetadata,OPᵢₙ<:AbstractOperator}
    isnothing(hermitian) && (hermitian = op_hermicity(in_operator))
    isnothing(float) && (float = op_float_type(in_operator))

    # Convert metadata
    in_metadata = op_metadata(in_operator)
    out_metadata = convert_metadata(
        Mₒᵤₜ,
        in_metadata;
        radii=radii,
        hermitian=hermitian,
        out_shconv=out_shconv,
        source_kwargs=source_kwargs,
        extra_kwargs=metadata_extra_kwargs,
    )

    # Build empty operator
    out_operator = build_operator(
        OPₒᵤₜ, out_metadata; value=zero(float), extra_kwargs=operator_extra_kwargs
    )

    # Convert operator data
    convert_data!(out_operator, in_operator)

    return out_operator
end

### DATA ###

"""
    convert_data!(out_operator, in_operator)

Transfer and transform data from `in_operator` into `out_operator` (mutating `out_operator`
in place). Currently dispatches on the `KeyedTrait` of both operators to select the correct
data types for conversion:

- `(NoKeydata, NoKeydata)` → route on `(Dₒᵤₜ, Dᵢₙ)`
- `(HasKeydata, NoKeydata)` → route on `(KDₒᵤₜ, Dₒᵤₜ, Dᵢₙ)`
- `(NoKeydata, HasKeydata)` → route on `(Dₒᵤₜ, KDᵢₙ, Dᵢₙ)`
- `(HasKeydata, HasKeydata)` → route on `(KDₒᵤₜ, Dₒᵤₜ, KDᵢₙ, Dᵢₙ)`

Concrete conversion methods are defined per format pair in `src/conversions/`.
"""
function convert_data!(
    out_operator::OPₒᵤₜ,
    in_operator::OPᵢₙ,
) where {
    OPₒᵤₜ<:AbstractOperator,
    OPᵢₙ<:AbstractOperator,
}
    return convert_data!(
        trait(KeyedTrait, OPₒᵤₜ),
        trait(KeyedTrait, OPᵢₙ),
        out_operator,
        in_operator,
    )
end

function convert_data!(
    ::NoKeydata,
    ::NoKeydata,
    out_operator::OPₒᵤₜ,
    in_operator::OPᵢₙ,
) where {
    OPₒᵤₜ<:AbstractOperator,
    OPᵢₙ<:AbstractOperator,
}
    Mₒᵤₜ = typeof(op_metadata(out_operator))
    Mᵢₙ = typeof(op_metadata(in_operator))
    Dₒᵤₜ = op_data_type(Mₒᵤₜ)
    Dᵢₙ = op_data_type(Mᵢₙ)
    return convert_data!(Dₒᵤₜ, Dᵢₙ, out_operator, in_operator)
end

function convert_data!(
    ::HasKeydata,
    ::NoKeydata,
    out_operator::OPₒᵤₜ,
    in_operator::OPᵢₙ,
) where {
    OPₒᵤₜ<:AbstractOperator,
    OPᵢₙ<:AbstractOperator,
}
    Mₒᵤₜ = typeof(op_metadata(out_operator))
    Mᵢₙ = typeof(op_metadata(in_operator))
    KDₒᵤₜ = op_keydata_type(Mₒᵤₜ)
    Dₒᵤₜ = op_data_type(Mₒᵤₜ)
    Dᵢₙ = op_data_type(Mᵢₙ)
    return convert_data!(KDₒᵤₜ, Dₒᵤₜ, Dᵢₙ, out_operator, in_operator)
end

function convert_data!(
    ::NoKeydata,
    ::HasKeydata,
    out_operator::OPₒᵤₜ,
    in_operator::OPᵢₙ,
) where {
    OPₒᵤₜ<:AbstractOperator,
    OPᵢₙ<:AbstractOperator,
}
    Mₒᵤₜ = typeof(op_metadata(out_operator))
    Mᵢₙ = typeof(op_metadata(in_operator))
    Dₒᵤₜ = op_data_type(Mₒᵤₜ)
    KDᵢₙ = op_keydata_type(Mᵢₙ)
    Dᵢₙ = op_data_type(Mᵢₙ)
    return convert_data!(Dₒᵤₜ, KDᵢₙ, Dᵢₙ, out_operator, in_operator)
end

function convert_data!(
    ::HasKeydata,
    ::HasKeydata,
    out_operator::OPₒᵤₜ,
    in_operator::OPᵢₙ,
) where {
    OPₒᵤₜ<:AbstractOperator,
    OPᵢₙ<:AbstractOperator,
}
    Mₒᵤₜ = typeof(op_metadata(out_operator))
    Mᵢₙ = typeof(op_metadata(in_operator))
    KDₒᵤₜ = op_keydata_type(Mₒᵤₜ)
    Dₒᵤₜ = op_data_type(Mₒᵤₜ)
    KDᵢₙ = op_keydata_type(Mᵢₙ)
    Dᵢₙ = op_data_type(Mᵢₙ)
    return convert_data!(KDₒᵤₜ, Dₒᵤₜ, KDᵢₙ, Dᵢₙ, out_operator, in_operator)
end

### METADATA ###

"""
    convert_metadata(::Type{Mₒᵤₜ}, in_metadata; kwargs...) -> AbstractMetadata

Convert metadata from one format to another. This is a three-stage pipeline:

1. [`convert_metadata_basic`](@ref) — convert source, SH convention, sparsity, basis.
2. `convert_metadata_extra` — convert extra fields (k-point, spins) via trait dispatch.
3. `convert_metadata_final` — assemble the concrete output metadata type.

`Mₒᵤₜ` may be a union type (e.g. `CanonicalBlockRealMetadata`); the concrete subtype is
chosen automatically based on which extra fields are present.

# Keyword arguments
- `radii=nothing`: build sparsity from neighbour list instead of converting.
- `hermitian=nothing`: override hermicity (defaults to input's).
- `out_shconv=nothing`: explicit output SH convention.
- `subbasis=nothing`: reduce basis to this subset of orbitals.
- `inverted=false`: if `true`, keep the complement of `subbasis`.
- `source_kwargs=NamedTuple()`: extra arguments for constructing the output source.
- `extra_kwargs=NamedTuple()`: extra arguments for trait-based conversion (e.g. `kpoint`).
"""
function convert_metadata(
    ::Type{Mₒᵤₜ}, in_metadata::AbstractMetadata;
    radii=nothing, hermitian=nothing, out_shconv=nothing,
    subbasis=nothing, inverted=false,
    source_kwargs=NamedTuple(), extra_kwargs=NamedTuple(),
) where {Mₒᵤₜ<:AbstractMetadata}

    # Construct basic metadata container
    out_basic_metadata = convert_metadata_basic(
        Mₒᵤₜ,
        in_metadata;
        radii=radii,
        hermitian=hermitian,
        out_shconv=out_shconv,
        source_kwargs=source_kwargs,
        subbasis=subbasis,
        inverted=inverted,
    )

    extra_args = convert_metadata_extra(
        Mₒᵤₜ,
        in_metadata,
        out_basic_metadata;
        subbasis=subbasis,
        inverted=inverted,
        extra_kwargs=extra_kwargs,
    )

    return convert_metadata_final(Mₒᵤₜ, out_basic_metadata, extra_args...)
end

### METADATA BASIC ###

"""
    convert_metadata_basic(::Type{Mₒᵤₜ}, in_metadata; kwargs...) -> BasicMetadataContainer

Convert the core metadata fields: derives the output source, computes the SH convention
delta, converts sparsity (optionally from radii), reorders the basis set to the output
SH convention, and optionally reduces it with a subbasis.
"""
function convert_metadata_basic(
    ::Type{Mₒᵤₜ}, in_metadata::Mᵢₙ;
    radii=nothing, hermitian=nothing, out_shconv=nothing,
    subbasis=nothing, inverted=false,
    source_kwargs=NamedTuple(),
) where {Mₒᵤₜ<:AbstractMetadata,Mᵢₙ<:AbstractMetadata}
    isnothing(hermitian) && (hermitian = op_hermicity(in_metadata))

    in_kind = op_kind(in_metadata)
    in_basisset = op_basisset(in_metadata)
    in_shconv = op_shconv(in_metadata)
    in_atoms = op_atoms(in_metadata)

    # Obtain output source
    Xₒᵤₜ = op_source_type(Mₒᵤₜ)
    out_source = Xₒᵤₜ(source_kwargs...)

    # Obtain spherical harmonics convention
    isnothing(out_shconv) && (out_shconv = default_shconv(out_source))
    Δshconv = out_shconv ∘ inv(in_shconv)

    # Convert sparsity
    Sₒᵤₜ = op_sparsity_type(Mₒᵤₜ)
    out_sparsity = convert_sparsity(Sₒᵤₜ, in_metadata; radii=radii, hermitian=hermitian)

    # Convert the spherical harmonics convention of the basis set
    out_basisset = convert_basisset_shconv(in_basisset, Δshconv)

    # Optionally reduce the basis with subbasis masks
    if !isnothing(subbasis)
        out_basisset = reduce_basisset(out_basisset, subbasis; inverted=inverted)
    end

    # Construct basic metadata container
    return BasicMetadataContainer(
        in_kind, out_source, out_sparsity, out_basisset, out_shconv, in_atoms
    )
end

"""
    convert_sparsity(::Type{Sₒᵤₜ}, in_metadata; radii=nothing, hermitian=nothing)

Metadata-level sparsity conversion wrapper. If `radii` is provided, builds a new sparsity
pattern from a neighbour list; otherwise converts the input metadata's existing sparsity to
type `Sₒᵤₜ`, optionally changing hermicity.
"""
function convert_sparsity(
    ::Type{Sₒᵤₜ}, in_metadata::AbstractMetadata;
    radii=nothing, hermitian=nothing,
) where {Sₒᵤₜ<:AbstractSparsity}
    isnothing(hermitian) && (hermitian = op_hermicity(in_metadata))
    if !isnothing(radii)
        in_atoms = op_atoms(in_metadata)
        return build_sparsity(Sₒᵤₜ, in_atoms, radii; hermitian=hermitian)
    end

    return convert_sparsity(Sₒᵤₜ, in_metadata, nothing; hermitian=hermitian)
end

function convert_sparsity(
    ::Type{Sₒᵤₜ}, in_metadata::AbstractMetadata, ::Nothing;
    hermitian=nothing,
) where {Sₒᵤₜ<:AbstractSparsity}
    isnothing(hermitian) && (hermitian = op_hermicity(in_metadata))
    in_sparsity = op_sparsity(in_metadata)
    in_basisset = op_basisset(in_metadata)

    return convert_sparsity(Sₒᵤₜ, in_sparsity, in_basisset; hermitian=hermitian)
end

### METADATA EXTRA ###

"""
    convert_metadata_extra(::Type{Mₒᵤₜ}, in_metadata, out_basic_metadata; kwargs...)

Convert extra metadata fields (k-point, spins) based on `extrafield_traittypes(Mₒᵤₜ)`.
Traits are sorted alphabetically and each dispatches to a trait-specific handler.
"""
function convert_metadata_extra(
    ::Type{Mₒᵤₜ}, in_metadata::AbstractMetadata, out_basic_metadata::BasicMetadataContainer;
    subbasis=nothing, inverted=false, extra_kwargs=NamedTuple(),
) where {Mₒᵤₜ<:AbstractMetadata}
    traits = extrafield_traittypes(Mₒᵤₜ)
    sorted_traits = view(traits, sortperm(Symbol.(traits)))
    sorted_extras = map(sorted_traits) do T
        convert_metadata_extra(
            T,
            Mₒᵤₜ,
            in_metadata,
            out_basic_metadata;
            subbasis=subbasis,
            inverted=inverted,
            extra_kwargs=extra_kwargs,
        )
    end
    return skipmissing(sorted_extras)
end

function convert_metadata_extra(
    ::Type{SpaceTrait}, ::Type{Mₒᵤₜ}, in_metadata::Mᵢₙ,
    out_basic_metadata::BasicMetadataContainer;
    subbasis=nothing, inverted=false, extra_kwargs=NamedTuple(),
) where {Mₒᵤₜ<:AbstractMetadata,Mᵢₙ<:AbstractMetadata}
    return convert_kpoint(
        Mₒᵤₜ, trait(SpaceTrait, Mᵢₙ), in_metadata; extra_kwargs=extra_kwargs
    )
end

function convert_kpoint(
    ::Type{Mₒᵤₜ}, ::RecipSpace, in_metadata::AbstractMetadata;
    extra_kwargs=NamedTuple(),
) where {Mₒᵤₜ<:AbstractMetadata}
    return op_kpoint(in_metadata)
end

function convert_kpoint(
    ::Type{Mₒᵤₜ}, ::RealSpace, in_metadata::AbstractMetadata;
    extra_kwargs=NamedTuple(),
) where {Mₒᵤₜ<:AbstractMetadata}
    @argcheck :kpoint ∈ keys(extra_kwargs)
    return extra_kwargs.kpoint
end

function convert_metadata_extra(
    ::Type{SpinTrait}, ::Type{Mₒᵤₜ}, in_metadata::Mᵢₙ,
    out_basic_metadata::BasicMetadataContainer;
    subbasis=nothing, inverted=false, extra_kwargs=NamedTuple(),
) where {Mₒᵤₜ<:AbstractMetadata,Mᵢₙ<:AbstractMetadata}
    return convert_spins(
        Mₒᵤₜ,
        trait(SpinTrait, Mᵢₙ),
        in_metadata,
        out_basic_metadata;
        subbasis=subbasis,
        inverted=inverted,
        extra_kwargs=extra_kwargs,
    )
end

# SpinTrait is attached to this union, but we know we need to construct RealMetadata in this case,
# which means:
# 1) Spins cannot be converted (nothing to convert)
# 2) Spins cannot be passed to extra_kwargs (no meaningful spins to attach)
# The way we alleviate this is by returning missing
function convert_metadata_extra(
    ::Type{SpinTrait}, ::Type{Mₒᵤₜ}, in_metadata::Mᵢₙ,
    out_basic_metadata::BasicMetadataContainer;
    subbasis=nothing, inverted=false, extra_kwargs=NamedTuple(),
) where {Mₒᵤₜ<:Union{<:RealMetadata,<:SpinRealMetadata},Mᵢₙ<:RealMetadata}
    return missing
end

# Same as above
function convert_metadata_extra(
    ::Type{SpinTrait}, ::Type{Mₒᵤₜ}, in_metadata::Mᵢₙ,
    out_basic_metadata::BasicMetadataContainer;
    subbasis=nothing, inverted=false, extra_kwargs=NamedTuple(),
) where {Mₒᵤₜ<:Union{<:RecipMetadata,<:SpinRecipMetadata},Mᵢₙ<:RecipMetadata}
    return missing
end

# Same as above
function convert_metadata_extra(
    ::Type{SpinTrait}, ::Type{Mₒᵤₜ}, in_metadata::Mᵢₙ,
    out_basic_metadata::BasicMetadataContainer;
    subbasis=nothing, inverted=false, extra_kwargs=NamedTuple(),
) where {Mₒᵤₜ<:Union{<:RecipMetadata,<:SpinRecipMetadata},Mᵢₙ<:RealMetadata}
    return missing
end

# Same as above
function convert_metadata_extra(
    ::Type{SpinTrait}, ::Type{Mₒᵤₜ}, in_metadata::Mᵢₙ,
    out_basic_metadata::BasicMetadataContainer;
    subbasis=nothing, inverted=false, extra_kwargs=NamedTuple(),
) where {Mₒᵤₜ<:Union{<:RealMetadata,<:SpinRealMetadata},Mᵢₙ<:RecipMetadata}
    return missing
end

function convert_spins(
    ::Type{Mₒᵤₜ}, ::NoSpin, in_metadata::AbstractMetadata,
    out_basic_metadata::BasicMetadataContainer;
    subbasis=nothing, inverted=false, extra_kwargs=NamedTuple(),
) where {Mₒᵤₜ<:AbstractMetadata}
    @argcheck :spins ∈ keys(extra_kwargs)
    return extra_kwargs.spins
end

function convert_spins(
    ::Type{Mₒᵤₜ}, ::HasSpin, in_metadata::AbstractMetadata,
    out_basic_metadata::BasicMetadataContainer;
    subbasis=nothing, inverted=false, extra_kwargs=NamedTuple(),
) where {Mₒᵤₜ<:AbstractMetadata}
    in_spins = op_spins(in_metadata)

    in_basisset = op_basisset(in_metadata)
    out_basisset = op_basisset(out_basic_metadata) # out_basisset is potentially a subbasisset already

    in_source = op_source(in_metadata)
    out_source = op_source(out_basic_metadata)

    in_shconv = op_shconv(in_metadata)
    out_shconv = op_shconv(out_basic_metadata)
    Δshconv = out_shconv ∘ inv(in_shconv)

    # Convert the spherical harmonics convention of spins
    if !isidentity(Δshconv)
        in_spins = convert_spins_shconv(in_spins, in_basisset, Δshconv)
    end

    # If subbasis is present, perform basis set reduction for spins.
    # The reduce_spins function accepts a basisset directly, but we compute masks here
    # instead because out_basisset would already be projected, meaning we cannot use it
    # to compute the masks. Furthermore, we cannot use in_basisset directly either
    # because it might not have the correct spherical harmonics convention
    if !isnothing(subbasis)
        in_basis = in_basisset.basis
        in_subbasis_masks = get_subbasis_masks(in_basisset, subbasis; inverted=inverted)
        out_subbasis_masks = convert_speciesdict_shconv(
            in_subbasis_masks, in_basis, Δshconv
        )
        in_spins = reduce_spins(in_spins, out_subbasis_masks)
    end

    # Make final changes due to the source change (e.g. reorder up and down spins if
    # the two sources don't agree). This often might leave the spins unchanged.
    if !isequal(in_source, out_source)
        in_spins = convert_spins_source(in_spins, out_basisset, in_source, out_source)
    end

    out_spins = in_spins
    return out_spins
end

### METADATA FINAL ###

"""
    convert_metadata_final(::Type{Mₒᵤₜ}, out_basic_metadata, extra_args...)

Assemble the final concrete metadata type from the basic metadata container and any extra
fields (k-point, spins). When `Mₒᵤₜ` is a union type, the concrete subtype is chosen based
on which extra arguments are present. For example:
- No extras → `RealMetadata`
- Spins only → `SpinRealMetadata`
- K-point only → `RecipMetadata`
- Both → `SpinRecipMetadata`
"""
function convert_metadata_final(
    ::Type{Mₒᵤₜ}, out_basic_metadata::BasicMetadataContainer
) where {Mₒᵤₜ<:Union{<:RealMetadata,<:SpinRealMetadata}}
    return RealMetadata(out_basic_metadata)
end

function convert_metadata_final(
    ::Type{Mₒᵤₜ}, out_basic_metadata::BasicMetadataContainer, out_spins::SpinsMetadata
) where {Mₒᵤₜ<:Union{<:RealMetadata,<:SpinRealMetadata}}
    return SpinRealMetadata(out_basic_metadata, out_spins)
end

function convert_metadata_final(
    ::Type{Mₒᵤₜ}, out_basic_metadata::BasicMetadataContainer, out_kpoint::SVector{3}
) where {Mₒᵤₜ<:Union{<:RecipMetadata,<:SpinRecipMetadata}}
    return RecipMetadata(out_basic_metadata, out_kpoint)
end

function convert_metadata_final(
    ::Type{Mₒᵤₜ}, out_basic_metadata::BasicMetadataContainer, out_kpoint::SVector{3},
    out_spins::SpinsMetadata
) where {Mₒᵤₜ<:Union{<:RecipMetadata,<:SpinRecipMetadata}}
    return SpinRecipMetadata(out_basic_metadata, out_spins, out_kpoint)
end

### ZERO OUT DATA ###

"""
    zero_out_data!(out_operator, in_metadata)

Zero out the data entries of `out_operator` whose sparsity is *not* contained in the sparsity
of `in_metadata`, mutating `out_operator` in place. Entries of `out_operator` that fall under
`in_metadata`'s sparsity keep their current values; every other entry is set to zero.

Unlike [`convert_operator`](@ref) / [`convert_data!`](@ref), which build a new operator with a
new (and therefore differently-typed) metadata, this routine changes only the data values and
leaves the metadata — and hence the operator's type — untouched. This is useful for restricting
an operator to a sparser support (e.g. a shorter-range neighbour pattern) without changing its
type, so that all methods defined for that type remain applicable.

Dispatches on the `KeyedTrait` of `out_operator`, then on the concrete output data type and the
input sparsity type. Concrete methods are defined per format in `src/conversions/`.

# Arguments
- `out_operator`: the operator whose data is zeroed in place. Its metadata is not modified.
- `in_metadata`: metadata whose sparsity determines which entries of `out_operator` are retained.
"""
function zero_out_data!(
    out_operator::OPₒᵤₜ,
    in_metadata::Mᵢₙ,
) where {
    OPₒᵤₜ<:AbstractOperator,
    Mᵢₙ<:AbstractMetadata,
}
    return zero_out_data!(
        trait(KeyedTrait, OPₒᵤₜ),
        out_operator,
        in_metadata,
    )
end

function zero_out_data!(
    ::NoKeydata,
    out_operator::OPₒᵤₜ,
    in_metadata::Mᵢₙ,
) where {
    OPₒᵤₜ<:AbstractOperator,
    Mᵢₙ<:AbstractMetadata,
}
    Mₒᵤₜ = typeof(op_metadata(out_operator))
    Dₒᵤₜ = op_data_type(Mₒᵤₜ)
    Sᵢₙ = op_sparsity_type(Mᵢₙ)
    return zero_out_data!(Dₒᵤₜ, Sᵢₙ, out_operator, in_metadata)
end

function zero_out_data!(
    ::HasKeydata,
    out_operator::OPₒᵤₜ,
    in_metadata::Mᵢₙ,
) where {
    OPₒᵤₜ<:AbstractOperator,
    Mᵢₙ<:AbstractMetadata,
}
    Mₒᵤₜ = typeof(op_metadata(out_operator))
    KDₒᵤₜ = op_keydata_type(Mₒᵤₜ)
    Dₒᵤₜ = op_data_type(Mₒᵤₜ)
    Sᵢₙ = op_sparsity_type(Mᵢₙ)
    return zero_out_data!(KDₒᵤₜ, Dₒᵤₜ, Sᵢₙ, out_operator, in_metadata)
end
