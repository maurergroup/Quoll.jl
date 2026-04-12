### FACTORIES ###

function Base.zero(operator::OP) where {OP<:AbstractOperator}
    return build_operator(
        OP,
        op_metadata(operator);
        value=zero(op_float_type(operator)),
        initialised=true,
    )
end

function Base.similar(operator::OP) where {OP<:AbstractOperator}
    return build_operator(
        OP,
        op_metadata(operator);
        type=op_float_type(operator),
        initialised=false,
    )
end

### OPERATOR ###

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
