### OPERATOR ###

# This function implies only metadata and data needs to be converted,
# everything else (e.g. keydata) can be built from metadata. Might not be true
# if operators have something like gradients attached
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

# Sparsity is normally used during conversion. This would imply that we need
# to use sparsity as trait too, but we can assume a data format implies a particular
# sparsity type.
# However, keydata does not imply sparsity, which is why we need data traits
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
    # to compute the masks. On the other hand, we cannot use in_basisset directly either
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
        # TODO: convert_spins_source doesn't exist currently
        in_spins = convert_spins_source(in_spins, out_basisset, in_source, out_source)
    end

    out_spins = in_spins
    return out_spins
end

### METADATA FINAL ###
# Often it's more convenient to use a union as the output metadata type
# because this doesn't require knowing whether the data has spin degrees of freedom
# or not during the call. Therefore, methods that accept unions were defined.
# However, those methods work for non-union Mₒᵤₜ types out of the box.

# Conversions for unions

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
    ::Type{Mₒᵤₜ}, out_basic_metadata::BasicMetadataContainer, out_spins::SpinsMetadata,
    out_kpoint::SVector{3}
) where {Mₒᵤₜ<:Union{<:RecipMetadata,<:SpinRecipMetadata}}
    return SpinRecipMetadata(out_basic_metadata, out_spins, out_kpoint)
end