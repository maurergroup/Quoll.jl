abstract type AbstractOperator end
abstract type AbstractOperatorMetadata end

function get_readformat end
function get_writeformat end

function get_avail_operatorkinds end
function get_avail_filenames end

# Find available operators in the supplied dir
function find_operatorkinds(dir::AbstractString, format::Type{<:AbstractOperator})
    found_operatorkinds = AbstractOperatorKind[]
    for operatorkind in get_avail_operatorkinds(format)
        names = get_avail_filenames(operatorkind, Val(operatorkind.tag), format)
        ps = joinpath.(dir, names)
        if any(ispath.(ps))
            push!(found_operatorkinds, operatorkind)
        end
    end
    return found_operatorkinds
end