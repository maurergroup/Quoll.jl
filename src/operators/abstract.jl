using Base

abstract type AbstractOperator end
abstract type AbstractOperatorMetadata end

function Base.show(io::IO, operator::AbstractOperator)
    print(io, nameof(typeof(operator)))
    print(io, "($(operator.kind))")
end

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