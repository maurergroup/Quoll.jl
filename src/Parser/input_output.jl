using Configurations
using ArgCheck
using ..Utils
using ..OperatorIO

@option struct InputParams <: AbstractQuollParams
    format::Type{<:AbstractOperator}
    directories::Vector{String}
    operators::Vector{AbstractOperatorKind} = [Hamiltonian(:ref), Hamiltonian(:pred), Overlap(:ref)]

    function InputParams(format, directories, operators)
        # Convert `dirs` to absolute paths
        directories = abspath.(directories)

        # Check if directories exist
        @argcheck all(isdir.(directories)) "Some of the supplied input directories do not exist"

        # Walk the directories to obtain all the paths
        directories = vcat(Utils.find_leafdirs.(directories)...)

        new(format, directories, operators)
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
