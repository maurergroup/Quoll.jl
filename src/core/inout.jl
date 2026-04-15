### LOADING ###

"""
    load_operators(::Type{OP}, ::Type{M}, dir, kinds) -> Vector{OP}

Load multiple operators from disk. Calls [`load_operator`](@ref) for each `OperatorKind` in
`kinds`.
"""
function load_operators(
    ::Type{OP}, ::Type{M}, dir::AbstractString, kinds::AbstractVector{<:OperatorKind}
) where {OP<:AbstractOperator,M<:AbstractMetadata}
    return load_operator.(OP, M, dir, kinds)
end

"""
    load_operator(::Type{OP}, ::Type{M}, dir, kind) -> OP

Load a single operator of type `OP` from directory `dir`. Loads metadata and data for the
given `OperatorKind`, then assembles via [`build_operator`](@ref). Format-specific loading
is dispatched through `load_metadata_basic`, `load_metadata`, and `load_data` methods
defined per format.
"""
function load_operator(
    ::Type{OP}, ::Type{M}, dir::AbstractString, kind::OperatorKind
) where {OP<:AbstractOperator,M<:AbstractMetadata}

    # Load metadata
    metadata = load_metadata(M, dir, kind)

    # Load data
    data = load_data(M, dir, kind)

    return build_operator(OP, metadata, data)
end

function load_metadata(
    ::Type{M}, dir::AbstractString, kind::OperatorKind;
) where {M<:AbstractMetadata}

    # Load basic metadata
    basic_metadata = load_metadata_basic(M, dir, kind)

    return load_metadata(M, dir, basic_metadata)
end

### WRITING ###

"""
    write_operators(::Type{M}, dir, operators)

Write a collection of operators to disk in the format determined by metadata type `M`.
Creates the directory if it doesn't exist.
"""
function write_operators end

### PARSING ###

"""
    find_operatorkinds(::Type{M}, dir) -> Vector{OperatorKind}

Scan directory `dir` for available operator data files matching metadata format `M`.
Returns the list of `OperatorKind`s whose files are found on disk.
"""
function find_operatorkinds(::Type{M}, dir::AbstractString) where {M<:AbstractMetadata}
    found_operatorkinds = OperatorKind[]
    for operatorkind in get_avail_operatorkinds(M)
        names = get_avail_filenames(M, operatorkind)
        ps = joinpath.(dir, names)
        if any(ispath.(ps))
            push!(found_operatorkinds, operatorkind)
        end
    end
    return found_operatorkinds
end

function get_avail_filenames(metadata::M) where {M<:AbstractMetadata}
    operatorkind = op_kind(metadata)
    return get_avail_filenames(M, operatorkind)
end

# Most of the time it's going to be only a single name.
# However, in some special cases we need to define specialised methods,
# e.g. in FHIaims we have both .h5 and .out files.
# Maybe it would be cleaner to have a single name; this could be achieved
# by having a field in FHIaimsFormat
function get_avail_filenames(
    ::Type{M}, operatorkind::OperatorKind
) where {M<:AbstractMetadata}
    return [get_avail_filename(M, statictuple(operatorkind))]
end