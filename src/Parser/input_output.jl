using Configurations
using ArgCheck
using AtomsBase
using Unitful
using UnitfulAtomic

using ..Utils
using ..OperatorIO

@option struct InputParams <: AbstractQuollParams
    format::Type{<:AbstractOperator}
    directories::Vector{String}
    operators::Vector{AbstractOperatorKind} = avail_operatorkinds(format)
    radii::Union{Dict{ChemicalSpecies, Quantity{Float64, Unitful.𝐋, typeof(u"Å")}}, Nothing} = nothing

    function InputParams(format, directories, operators, radii)
        # Convert `dirs` to absolute paths
        directories = abspath.(directories)

        # Check if directories exist
        @argcheck all(isdir.(directories)) "Some of the supplied input directories do not exist"

        # Walk the directories to obtain all the paths
        directories = vcat(Utils.find_leafdirs.(directories)...)

        new(format, directories, operators, radii)
    end
end

@option struct OutputParams <: AbstractQuollParams
    format::Type{<:AbstractOperator}
    directory::String
    
    function OutputParams(format, directory)
        directory = abspath(directory)
        new(format, directory)
    end
end

function Configurations.from_dict(
    ::Type{InputParams},
    ::OptionField{:format},
    ::Type{Type{<:AbstractOperator}},
    s,
)
    @argcheck isa(s, String)
    
    symbol = Symbol(Utils.normalize_comparison(s))
    @argcheck hasmethod(read_format, Tuple{Val{symbol}}) "Writing matrix in format $s is unavailable or the format doesn't exist"

    return read_format(Val(symbol))
end

function Configurations.from_dict(
    ::Type{OutputParams},
    ::OptionField{:format},
    ::Type{Type{<:AbstractOperator}},
    s,
)
    @argcheck isa(s, String)
    
    symbol = Symbol(Utils.normalize_comparison(s))
    @argcheck hasmethod(write_format, Tuple{Val{symbol}}) "Writing matrix in format $s is unavailable or the format doesn't exist"

    return write_format(Val(symbol))
end

function Configurations.from_dict(
    ::Type{InputParams},
    ::OptionField{:directories},
    ::Type{Vector{String}},
    dirs,
)
    # Ensure `dirs` is of type Vector{String}
    @argcheck isa(dirs, String) || isa(dirs, Vector{String})
    return isa(dirs, String) ? [dirs] : dirs
end

function Configurations.from_dict(
    ::Type{InputParams},
    ::OptionField{:operators},
    ::Type{Vector{AbstractOperatorKind}},
    operators
)
    # Ensure `operators` is of type Vector{String}
    @argcheck isa(operators, String) || isa(operators, Vector{String})
    isa(operators, String) && (operators = [operators])

    valid_symbols = []
    for operator in operators
        symbol = Symbol(Utils.normalize_comparison(operator))
        @argcheck hasmethod(operator_kind, Tuple{Val{symbol}}) "Operator kind $operator is unavailable or doesn't exist"
        push!(valid_symbols, symbol)
    end

    return operator_kind.(Val.(valid_symbols))
end

function parse_radius(radius_string)
    re = r"\s*(\w+)\s*\+?(\d*\.?\d+)\s*(\w+)?\s*"
    radius_match = match(re, radius_string)

    z = ChemicalSpecies(Symbol(radius_match.captures[1]))
    r = parse(Float64, radius_match.captures[2])
    unit = radius_match.captures[3]
    
    if unit === nothing
        @warn "Units to radii field not supplied, assuming Å"
        r = r * u"Å"
    else
        unit = Utils.normalize_comparison(unit)
        if unit == "bohr"
            r = uconvert(u"Å", r * u"bohr")
        elseif unit == "ang"
            r = r * u"Å"
        else
            throw(ArgumentError("Unknown radius unit $unit"))
        end
    end
    
    return z, r
end

function Configurations.from_dict(
    ::Type{InputParams},
    ::OptionField{:radii},
    ::Type{Dict{ChemicalSpecies, Quantity{Float64, Unitful.𝐋, typeof(u"Å")}}},
    radii
)
    # Ensure `radii` is of type Vector{String}
    @argcheck isa(radii, String) || isa(radii, Vector{String})
    isa(radii, String) && (radii = [radii])

    return Dict(Pair(parse_radius(radius)...) for radius in radii)
end
