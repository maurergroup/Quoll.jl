using Base

abstract type AbstractOperator end
abstract type AbstractOperatorMetadata end

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

get_atoms(metadata::AbstractOperatorMetadata) = metadata.atoms
get_sparsity(metadata::AbstractOperatorMetadata) = metadata.sparsity
get_basisset(metadata::AbstractOperatorMetadata) = metadata.basisset
get_spins(metadata::AbstractOperatorMetadata) = metadata.spins

# TODO: Even though the fields of metadata might be the same among different formats,
# the constructors will differ.
# An alternative would be to use dispatch on the format instead,
# e.g. OperatorMetadata(..., ::Type{<:AbstractOperator}).
# For now I will keep metadata as separate structs for each format,
# but if I find that fields are the same I should consider refactoring
# the code in favour of dispatch on operator format (more consistent
# with the rest and fewer unnecessary types).
# (N.B. in that case OperatorMetadata would become a more parametric type
# than before because we wouldn't know the sparsity type)

# Again things related to parsing below

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
    pairtypes = (Pair(Val(pair[1]), Val(pair[2])) for pair in pairs(operatorkind.tags))
    return get_avail_filenames(operatorkind, pairtypes..., T)
end

# Could implement basic versions of those like for convert_operator
function load_operators end
function write_operators end

# Could be specialised for particular operator conversions if the conversion can be done in a more efficient way
# e.g. if metadata is always shared between operators
function convert_operators(in_operators, out_operator_type::Type{<:AbstractOperator};
    radii = nothing, hermitian = nothing, float = nothing)
    return convert_operator.(in_operators, out_operator_type; radii = radii, hermitian = hermitian, float = float)
end

# Requires defining appropriate constructor for T2.
# This function is just an alias for the constructor
function convert_operator(in_operator::T1, ::Type{T2};
    radii = nothing, hermitian = nothing, float = nothing) where {T1<:AbstractOperator, T2<:AbstractOperator}

    isnothing(hermitian) && (hermitian = get_sparsity(in_operator).hermitian)
    isnothing(float) && (float = get_float(in_operator))

    return T2(in_operator, radii = radii, hermitian = hermitian, float = float)
end