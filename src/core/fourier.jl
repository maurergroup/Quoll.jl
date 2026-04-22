### FOURIER TRANSFORM ###

"""
    fourier_transform(::Type{OPₒᵤₜ}, ::Type{Mₒᵤₜ}, in_operator, kpoints; kwargs...)

Fourier-transform a real-space operator to reciprocal space at one or more k-points.
Returns a vector of reciprocal-space operators (one per k-point).

Internally: precomputes phases from the input sparsity's image vectors, converts metadata
(injecting the k-point), allocates a complex output operator, and calls
[`fourier_transform_data!`](@ref) for each k-point.

# Arguments
- `OPₒᵤₜ`: output operator type (`<:AbstractOperator`).
- `Mₒᵤₜ`: output metadata type (e.g. `CanonicalDenseRecipMetadata`).
- `in_operator`: real-space input operator.
- `kpoints`: a single k-point vector or a collection of k-point vectors.

# Keyword arguments
- `out_shconv=nothing`: explicit SH convention for the output (defaults to the source's).
- `source_kwargs=NamedTuple()`: extra keyword arguments for constructing the output source.
"""
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

    # Add the k-point to the extra metadata keyword arguments. The k-point is also
    # converted to a static array (this is what `convert_metadata_final()` expects;
    # if `kpoint` is already static, this doesn't have an effect)
    kpoint_static = SVector{3}(kpoint)
    ext_metadata_extra_kwargs = merge((; kpoint=kpoint_static), metadata_extra_kwargs)

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

"""
    fourier_transform_data!(out_operator, in_operator, phases_k)

Perform the forward Fourier transform of data from `in_operator` into `out_operator`
(mutating in place) using precomputed phase factors `phases_k`. Currently dispatches on
`KeyedTrait` of both operators, mirroring the 4-way dispatch pattern of [`convert_data!`](@ref).

Concrete implementations are in `src/conversions/`.
"""
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

"""
    inv_fourier_transform_data!(out_operator, in_operator, phases_k, weight)

Perform the inverse Fourier transform, accumulating the contribution of `in_operator`
(at a single k-point) into the real-space `out_operator`, weighted by `weight`.
Designed to be called in a loop over k-points, summing weighted contributions.

Dispatches on `KeyedTrait` of both operators, same pattern as `fourier_transform_data!`.
"""
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