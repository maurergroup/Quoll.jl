### FOURIER TRANSFORM ###

function fourier_transform(
    ::Type{OPₒᵤₜ}, ::Type{Mₒᵤₜ}, in_operator::AbstractOperator, kpoints;
    out_shconv=nothing, source_kwargs=NamedTuple(),
) where {OPₒᵤₜ<:AbstractOperator,Mₒᵤₜ<:AbstractMetadata}
    if !isa(kpoints, AbstractVector{<:AbstractVector})
        kpoints = [kpoints]
    end
    in_images = op_images(op_sparsity(in_operator))
    phases = eachcol(precompute_phases(kpoints, in_images))
    return fourier_transform.(
        OPₒᵤₜ, Mₒᵤₜ, Ref(in_operator), kpoints, phases;
        out_shconv=out_shconv, source_kwargs=source_kwargs,
    )
end

function fourier_transform(
    ::Type{OPₒᵤₜ}, ::Type{Mₒᵤₜ}, in_operator::AbstractOperator, kpoint, phases_k;
    out_shconv=nothing, source_kwargs=NamedTuple(),
    operator_extra_kwargs=NamedTuple(), metadata_extra_kwargs=NamedTuple(),
) where {OPₒᵤₜ<:AbstractOperator,Mₒᵤₜ<:AbstractMetadata}
    hermitian = op_hermicity(in_operator)

    # Add the k-point to the extra metadata keyword arguments
    ext_metadata_extra_kwargs = merge((; kpoint=kpoint), metadata_extra_kwargs)

    # Convert metadata
    in_metadata = op_metadata(in_operator)
    out_metadata = convert_metadata(
        Mₒᵤₜ,
        in_metadata;
        hermitian=hermitian,
        out_shconv=out_shconv,
        source_kwargs=source_kwargs,
        extra_kwargs=ext_metadata_extra_kwargs,
    )

    return fourier_transform(
        OPₒᵤₜ, in_operator, out_metadata, phases_k; extra_kwargs=operator_extra_kwargs
    )
end

function fourier_transform(
    ::Type{OPₒᵤₜ}, in_operator::AbstractOperator, out_metadata::AbstractMetadata, phases_k;
    extra_kwargs=NamedTuple(),
) where {OPₒᵤₜ<:AbstractOperator}
    float = op_float_type(in_operator)
    cfloat = float <: Complex ? float : Complex{float}

    # Build empty operator
    out_operator = build_operator(
        OPₒᵤₜ, out_metadata; type=cfloat, extra_kwargs=extra_kwargs
    )

    # Perform Fourier transform
    fourier_transform_data!(out_operator, in_operator, phases_k)

    return out_operator
end

function fourier_transform_data!(
    out_operator::OPₒᵤₜ,
    in_operator::OPᵢₙ,
    phases_k,
) where {
    OPₒᵤₜ<:AbstractOperator,
    OPᵢₙ<:AbstractOperator,
}
    return fourier_transform_data!(
        trait(KeyedTrait, OPₒᵤₜ),
        trait(KeyedTrait, OPᵢₙ),
        out_operator,
        in_operator,
        phases_k,
    )
end

function fourier_transform_data!(
    ::NoKeydata,
    ::NoKeydata,
    out_operator::OPₒᵤₜ,
    in_operator::OPᵢₙ,
    phases_k,
) where {
    OPₒᵤₜ<:AbstractOperator,
    OPᵢₙ<:AbstractOperator,
}
    Mₒᵤₜ = typeof(op_metadata(out_operator))
    Mᵢₙ = typeof(op_metadata(in_operator))
    Dₒᵤₜ = op_data_type(Mₒᵤₜ)
    Dᵢₙ = op_data_type(Mᵢₙ)
    return fourier_transform_data!(Dₒᵤₜ, Dᵢₙ, out_operator, in_operator, phases_k)
end

function fourier_transform_data!(
    ::HasKeydata,
    ::NoKeydata,
    out_operator::OPₒᵤₜ,
    in_operator::OPᵢₙ,
    phases_k,
) where {
    OPₒᵤₜ<:AbstractOperator,
    OPᵢₙ<:AbstractOperator,
}
    Mₒᵤₜ = typeof(op_metadata(out_operator))
    Mᵢₙ = typeof(op_metadata(in_operator))
    KDₒᵤₜ = op_keydata_type(Mₒᵤₜ)
    Dₒᵤₜ = op_data_type(Mₒᵤₜ)
    Dᵢₙ = op_data_type(Mᵢₙ)
    return fourier_transform_data!(KDₒᵤₜ, Dₒᵤₜ, Dᵢₙ, out_operator, in_operator, phases_k)
end

function fourier_transform_data!(
    ::NoKeydata,
    ::HasKeydata,
    out_operator::OPₒᵤₜ,
    in_operator::OPᵢₙ,
    phases_k,
) where {
    OPₒᵤₜ<:AbstractOperator,
    OPᵢₙ<:AbstractOperator,
}
    Mₒᵤₜ = typeof(op_metadata(out_operator))
    Mᵢₙ = typeof(op_metadata(in_operator))
    Dₒᵤₜ = op_data_type(Mₒᵤₜ)
    KDᵢₙ = op_keydata_type(Mᵢₙ)
    Dᵢₙ = op_data_type(Mᵢₙ)
    return fourier_transform_data!(Dₒᵤₜ, KDᵢₙ, Dᵢₙ, out_operator, in_operator, phases_k)
end

function fourier_transform_data!(
    ::HasKeydata,
    ::HasKeydata,
    out_operator::OPₒᵤₜ,
    in_operator::OPᵢₙ,
    phases_k,
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
    return fourier_transform_data!(
        KDₒᵤₜ, Dₒᵤₜ, KDᵢₙ, Dᵢₙ, out_operator, in_operator, phases_k
    )
end

### INVERSE FOURIER TRANSFORM ###

function inv_fourier_transform_data!(
    out_operator::OPₒᵤₜ,
    in_operator::OPᵢₙ,
    phases_k,
    weight,
) where {
    OPₒᵤₜ<:AbstractOperator,
    OPᵢₙ<:AbstractOperator,
}
    return inv_fourier_transform_data!(
        trait(KeyedTrait, OPₒᵤₜ),
        trait(KeyedTrait, OPᵢₙ),
        out_operator,
        in_operator,
        phases_k,
        weight,
    )
end

function inv_fourier_transform_data!(
    ::NoKeydata,
    ::NoKeydata,
    out_operator::OPₒᵤₜ,
    in_operator::OPᵢₙ,
    phases_k,
    weight,
) where {
    OPₒᵤₜ<:AbstractOperator,
    OPᵢₙ<:AbstractOperator,
}
    Mₒᵤₜ = typeof(op_metadata(out_operator))
    Mᵢₙ = typeof(op_metadata(in_operator))
    Dₒᵤₜ = op_data_type(Mₒᵤₜ)
    Dᵢₙ = op_data_type(Mᵢₙ)
    return inv_fourier_transform_data!(
        Dₒᵤₜ, Dᵢₙ, out_operator, in_operator, phases_k, weight
    )
end

function inv_fourier_transform_data!(
    ::HasKeydata,
    ::NoKeydata,
    out_operator::OPₒᵤₜ,
    in_operator::OPᵢₙ,
    phases_k,
    weight,
) where {
    OPₒᵤₜ<:AbstractOperator,
    OPᵢₙ<:AbstractOperator,
}
    Mₒᵤₜ = typeof(op_metadata(out_operator))
    Mᵢₙ = typeof(op_metadata(in_operator))
    KDₒᵤₜ = op_keydata_type(Mₒᵤₜ)
    Dₒᵤₜ = op_data_type(Mₒᵤₜ)
    Dᵢₙ = op_data_type(Mᵢₙ)
    return inv_fourier_transform_data!(
        KDₒᵤₜ, Dₒᵤₜ, Dᵢₙ, out_operator, in_operator, phases_k, weight
    )
end

function inv_fourier_transform_data!(
    ::NoKeydata,
    ::HasKeydata,
    out_operator::OPₒᵤₜ,
    in_operator::OPᵢₙ,
    phases_k,
    weight,
) where {
    OPₒᵤₜ<:AbstractOperator,
    OPᵢₙ<:AbstractOperator,
}
    Mₒᵤₜ = typeof(op_metadata(out_operator))
    Mᵢₙ = typeof(op_metadata(in_operator))
    Dₒᵤₜ = op_data_type(Mₒᵤₜ)
    KDᵢₙ = op_keydata_type(Mᵢₙ)
    Dᵢₙ = op_data_type(Mᵢₙ)
    return inv_fourier_transform_data!(
        Dₒᵤₜ, KDᵢₙ, Dᵢₙ, out_operator, in_operator, phases_k, weight
    )
end

function inv_fourier_transform_data!(
    ::HasKeydata,
    ::HasKeydata,
    out_operator::OPₒᵤₜ,
    in_operator::OPᵢₙ,
    phases_k,
    weight,
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
    return inv_fourier_transform_data!(
        KDₒᵤₜ, Dₒᵤₜ, KDᵢₙ, Dᵢₙ, out_operator, in_operator, phases_k, weight
    )
end