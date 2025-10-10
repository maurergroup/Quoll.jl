abstract type AbstractOperator end
abstract type AbstractOperatorMetadata end

function read_format end
function write_format end

function avail_operatorkinds end
function avail_filenames end

# Find available operators in the supplied dir
function find_operatorkinds(dir::AbstractString, format::Type{<:AbstractOperator})
    found_operatorkinds = AbstractOperatorKind[]
    for operatorkind in avail_operatorkinds(format)
        names = avail_filenames(operatorkind, Val(operatorkind.tag), format)
        ps = joinpath.(dir, names)
        if any(ispath.(ps))
            push!(found_operatorkinds, operatorkind)
        end
    end
    return found_operatorkinds
end