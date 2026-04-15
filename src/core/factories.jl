### FACTORIES ###

"""
    zero(operator::AbstractOperator)

Return a new operator of the same type and metadata, with all data elements set to zero.
"""
function Base.zero(operator::OP) where {OP<:AbstractOperator}
    return build_operator(
        OP,
        op_metadata(operator);
        value=zero(op_float_type(operator)),
        initialised=true,
    )
end

"""
    similar(operator::AbstractOperator)

Return a new operator of the same type and metadata, with allocated but uninitialised data.
"""
function Base.similar(operator::OP) where {OP<:AbstractOperator}
    return build_operator(
        OP,
        op_metadata(operator);
        type=op_float_type(operator),
        initialised=false,
    )
end

### OPERATOR ###

"""
    build_operator(::Type{OP}, metadata; kwargs...) -> OP
    build_operator(::Type{OP}, metadata, data; extra_kwargs=NamedTuple()) -> OP

Construct an operator of type `OP` (`Operator` or `KeyedOperator`) from `metadata`.

The first form allocates and initialises data according to the metadata's format. The second
form uses pre-existing `data` directly.

For `KeyedOperator`, keydata is automatically built from the data via `build_keydata`.

# Keyword arguments
- `value=0.0`: fill value for data arrays (determines element type if `type` is not given).
- `type=Nothing`: explicit element type override (e.g. `Float32`, `ComplexF64`).
- `initialised=true`: if `false`, allocate data without initialising (requires `type`).
- `subbasis=nothing`: if provided, reduce the basis set to this subbasis before building.
- `inverted=false`: if `true`, keep the complement of `subbasis` instead.
- `extra_kwargs=NamedTuple()`: passed through to `build_operator_extra`.
"""
function build_operator(
    ::Type{OP}, metadata::M;
    value=0.0, type=Nothing, initialised=true, subbasis=nothing, inverted=false,
    extra_kwargs=NamedTuple(),
) where {OP<:AbstractOperator,M<:AbstractMetadata}
    if !isnothing(subbasis)
        metadata = convert_metadata(M, metadata; subbasis=subbasis, inverted=inverted)
    end

    # Build data
    data = build_data(metadata; value=value, type=type, initialised=initialised)

    return build_operator(OP, metadata, data; extra_kwargs=extra_kwargs)
end

function build_operator(
    ::Type{OP}, metadata::AbstractMetadata, data::DataContainer;
    extra_kwargs=NamedTuple(),
) where {OP<:AbstractOperator}

    # Build extra args
    extra_args = build_operator_extra(OP, metadata, data; extra_kwargs=extra_kwargs)

    return build_operator_final(OP, metadata, data, extra_args...)
end

### OPERATOR EXTRA ###

"""
    build_operator_extra(::Type{OP}, metadata, data; extra_kwargs) -> extras

Build extra operator fields (e.g. keydata) required by operator type `OP`, determined by
its `extrafield_traittypes`. Returns a `skipmissing` iterator over the built extras.
"""
function build_operator_extra(
    ::Type{OP}, metadata::AbstractMetadata, data::DataContainer; extra_kwargs=NamedTuple()
) where {OP<:AbstractOperator}
    traits = extrafield_traittypes(OP)
    sorted_traits = view(traits, sortperm(Symbol.(traits)))
    sorted_extras = map(sorted_traits) do T
        build_operator_extra(T, metadata, data; extra_kwargs=extra_kwargs)
    end
    return skipmissing(sorted_extras)
end

function build_operator_extra(
    ::Type{KeyedTrait}, metadata::M, data::DataContainer;
    extra_kwargs=NamedTuple(),
) where {M<:AbstractMetadata}
    KD = op_keydata_type(M)
    return build_keydata(KD, metadata, data)
end

### OPERATOR FINAL ###

function build_operator_final(
    ::Type{<:Operator},
    metadata::AbstractMetadata,
    data::DataContainer,
)
    return Operator(metadata, data)
end

function build_operator_final(
    ::Type{<:KeyedOperator},
    metadata::AbstractMetadata,
    data::DataContainer,
    keydata::DataContainer,
)
    return KeyedOperator(metadata, data, keydata)
end

### DATA ###

"""
    build_data(metadata; value=0.0, type=Nothing, initialised=true)

Allocate a [`DataContainer`](@ref) matching the format described by `metadata`. Dispatches
to format-specific `build_data` methods based on `op_data_type(M)`.
"""
function build_data(
    metadata::M;
    value=0.0, type=Nothing, initialised=true,
) where {M<:AbstractMetadata}

    # Synchronise types of `value` and `type`
    if !(type <: Nothing)
        value = convert(type, value)
    else
        type = typeof(value)
    end

    # Type of uninitialised data needs to be known
    !initialised && @argcheck !(type <: Nothing)

    return build_data(op_data_type(M), metadata, value, initialised)
end
