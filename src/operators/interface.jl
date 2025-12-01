### METADATA INTERFACE ###

abstract type AbstractOperatorMetadata{A, S, B, P} end

# Accessors for properties which should be present in the operator metadata
get_atoms(metadata::AbstractOperatorMetadata) = metadata.atoms
get_sparsity(metadata::AbstractOperatorMetadata) = metadata.sparsity
get_basisset(metadata::AbstractOperatorMetadata) = metadata.basisset
get_spins(metadata::AbstractOperatorMetadata) = metadata.spins

### OPERATOR INTERFACE ###

abstract type AbstractOperator{O, T, D, M} end

function Base.show(io::IO, operator::AbstractOperator)
    print(io, nameof(typeof(operator)))
    print(io, "($(operator.kind))")
end

# Accessors for properties which should be present in the operator
get_kind(operator::AbstractOperator) = operator.kind
get_data(operator::AbstractOperator) = operator.data
get_metadata(operator::AbstractOperator) = operator.metadata

get_atoms(operator::AbstractOperator) = operator.metadata.atoms
get_sparsity(operator::AbstractOperator) = operator.metadata.sparsity
get_basisset(operator::AbstractOperator) = operator.metadata.basisset
get_spins(operator::AbstractOperator) = operator.metadata.spins

get_float(::AbstractOperator{O, T}) where {O, T} = T

### FACTORIES ###

function Base.zero(operator::T) where T<:AbstractOperator
    return build_operator(
        T, get_kind(operator), get_metadata(operator),
        uninit = false,
        value = zero(get_float(operator)),
    )
end

function Base.similar(operator::T) where T<:AbstractOperator
    return build_operator(
        T, get_kind(operator), get_metadata(operator),
        uninit = true,
        type = get_float(operator)
    )
end

function build_operator(::Type{T}, operatorkind::OperatorKind, metadata::AbstractOperatorMetadata;
    uninit = false, value = 0.0, type = Nothing) where T<:AbstractOperator
    constructor = Base.typename(T).wrapper
    return constructor(operatorkind, metadata, uninit = uninit, value = value, type = type)
end

function build_operator(::Type{T}, operatorkind::OperatorKind, data, metadata::AbstractOperatorMetadata
    ) where T<:AbstractOperator
    constructor = Base.typename(T).wrapper
    return constructor(operatorkind, data, metadata)
end

### BASIS SET PROJECTION ###

# Assumes only basis and spins explicitly depend on specific orbitals.
# If this is not the case for a given metadata type then this method
# needs to be specialised for that typ
function construct_subbasis_metadata(metadata::AbstractOperatorMetadata, subbasis::Vector{T};
    inverted = false) where T<:BasisMetadata

    basisset = get_basisset(metadata)
    spins = get_spins(metadata)

    subbasis_masks = get_subbasis_masks(basisset, subbasis, inverted = inverted)
    subbasisset = BasisSetMetadata(basisset, subbasis_masks)
    isnothing(spins) ? subspins = nothing : subspins = SpinsMetadata(spins, subbasis_masks)

    return construct_subbasis_metadata(metadata, subbasisset, subspins)
end

function construct_subbasis_metadata(metadata::M, subbasisset::BasisSetMetadata, subspins::Union{SpinsMetadata, Nothing}
    ) where M<:AbstractOperatorMetadata

    field_symbols = fieldnames(M)
    nfields = length(field_symbols)

    newfields = ntuple(Val(nfields)) do i_field
        field_symbol = field_symbols[i_field]
        if isequal(field_symbol, :basisset)
            return subbasisset
        elseif isequal(field_symbol, :spins)
            return subspins
        else
            getfield(metadata, field_symbol)
        end
    end

    return M(newfields...)
end

### WRITE INTERFACE ###

function write_operators(dir::AbstractString, operators)
    ispath(dir) || mkpath(dir)
    return write_operator.(dir, operators)
end

# Here using a constructor as below would not make sense
function write_operator end

### LOAD INTERFACE ###

function load_operators(dir::AbstractString, operatorkinds, operator_format::Type{<:AbstractOperator})
    return load_operator.(dir, operatorkinds, operator_format)
end

# TODO: if this is the API I am going for then I should define an appropriate constructor for FHIaims
# even though I defined a custom `load_operators`
function load_operator(dir::AbstractString, operatorkind::OperatorKind, ::Type{T}) where T<:AbstractOperator
    return T(dir, operatorkind)
end

### SPIN INTERFACE ###

# Constructors for SpinsMetadata from OperatorKind, BasisSetMetadata, and AbstractOperator.
# AbstractOperator is required for the SOC case. This is because for now in the case of SOC
# BasisSetMetadata.basis stores 2x the number of basis functions, and we don't know which of them
# are spin up and spin down. This type of information might depend on the format

function SpinsMetadata(kind::OperatorKind, basisset::BasisSetMetadata, format::Type{<:AbstractOperator})
    return SpinsMetadata(kind, basisset.basis, format)
end

function SpinsMetadata(kind::Hamiltonian, basis::SpeciesAnyDict, format::Type{<:AbstractOperator})
    return SpinsMetadata(kind, Val(kind.spin), basis, format)
end

function SpinsMetadata(::Hamiltonian, ::Val{:up}, basis::SpeciesAnyDict, ::Type{<:AbstractOperator})
    spins_dict = Base.ImmutableDict((
        z => fill(⬆, length(atom_basis))
        for (z, atom_basis) in pairs(basis)
    )...)
    return SpinsMetadata(spins_dict, false)
end

function SpinsMetadata(::Hamiltonian, ::Val{:down}, basis::SpeciesAnyDict, ::Type{<:AbstractOperator})
    spins_dict = Base.ImmutableDict((
        z => fill(⬇, length(atom_basis))
        for (z, atom_basis) in pairs(basis)
    )...)
    return SpinsMetadata(spins_dict, false)
end

SpinsMetadata(::Hamiltonian, ::Val{:none}, ::SpeciesAnyDict, ::Type{<:AbstractOperator}) = nothing

SpinsMetadata(::Overlap, ::SpeciesAnyDict, ::Type{<:AbstractOperator}) = nothing

### PARSING INTERFACE ###

function get_readformat end
function get_writeformat end

function get_avail_operatorkinds end
function get_avail_filenames end

# Find available operators in the supplied dir
function find_operatorkinds(dir::AbstractString, format::Type{<:AbstractOperator})
    found_operatorkinds = OperatorKind[]
    for operatorkind in get_avail_operatorkinds(format)
        names = get_avail_filenames(operatorkind, format)
        ps = joinpath.(dir, names)
        if any(ispath.(ps))
            push!(found_operatorkinds, operatorkind)
        end
    end
    return found_operatorkinds
end

function get_avail_filenames(operatorkind::OperatorKind, T::Type{<:AbstractOperator})
    pairtypes = Tuple(Pair(Val(pair[1]), Val(pair[2])) for pair in pairs(operatorkind.tags))
    return get_avail_filenames(operatorkind, pairtypes, T)
end

function get_avail_filenames(operatorkind::OperatorKind, pairtypes::Tuple, T::Type{<:AbstractOperator})
    return [get_avail_filename(operatorkind, pairtypes..., T)]
end

function get_avail_filename(operator::T) where T<:AbstractOperator
    operatorkind = get_kind(operator)
    pairtypes = (Pair(Val(pair[1]), Val(pair[2])) for pair in pairs(operatorkind.tags))
    return get_avail_filename(operatorkind, pairtypes..., T)
end