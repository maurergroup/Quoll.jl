using Configurations
using ArgCheck
using ..Utils
using ..OperatorIO

@option struct InputParams <: AbstractQuollParams
    format::Type{<:AtomBasisOperator}
    directories::Vector{String}

    function InputParams(format, directories)
        # Convert `dirs` to absolute paths
        directories = abspath.(directories)

        # Check if directories exist
        @argcheck all(isdir.(directories)) "Some of the supplied input directories do not exist"

        # Walk the directories to obtain all the paths
        directories = vcat(Utils.find_leafdirs.(directories)...)

        new(format, directories)
    end
end

@option struct OutputParams <: AbstractQuollParams
    format::Type{<:AtomBasisOperator}
    directory::String
    
    function OutputParams(format, directory)
        directory = abspath(directory)
        new(format, directory)
    end
end

function Configurations.from_dict(
    ::Type{InputParams},
    ::OptionField{:format},
    ::Type{Type{<:AtomBasisOperator}},
    s,
)
    @argcheck isa(s, String)
    
    symbol = Symbol(Utils.normalize_comparison(s))
    @argcheck hasmethod(implemented_read_format, Tuple{Val{symbol}}) "Writing matrix in format $s is unavailable or the format doesn't exist"

    return implemented_read_format(Val(symbol))
end

function Configurations.from_dict(
    ::Type{OutputParams},
    ::OptionField{:format},
    ::Type{Type{<:AtomBasisOperator}},
    s,
)
    @argcheck isa(s, String)
    
    symbol = Symbol(Utils.normalize_comparison(s))
    @argcheck hasmethod(implemented_write_format, Tuple{Val{symbol}}) "Writing matrix in format $s is unavailable or the format doesn't exist"

    return implemented_write_format(Val(symbol))
end

# TODO: Remove if we don't use REGISTRY approach
# function Configurations.from_dict(
#     ::Type{InputParams},
#     ::OptionField{:format},
#     ::Type{Type{<:AtomBasisOperator}},
#     format,
# )
#     @argcheck isa(format, String)

#     itype = findfirst(
#         x -> x == Utils.normalize_comparison(format),
#         lowercase.(string.(nameof.(READ_REGISTRY)))
#     )
#     @argcheck itype !== nothing "Reading matrix in format $format is unavailable or it doesn't exist"

#     return READ_REGISTRY[itype]
# end
#
# function Configurations.from_dict(
#     ::Type{OutputParams},
#     ::OptionField{:format},
#     ::Type{Type{<:AtomBasisOperator}},
#     format,
# )
#     @argcheck isa(format, String)

#     itype = findfirst(
#         x -> x == Utils.normalize_comparison(format),
#         lowercase.(string.(nameof.(WRITE_REGISTRY)))
#     )
#     @argcheck itype !== nothing "Writing matrix in format $format is unavailable or it doesn't exist"

#     return WRITE_REGISTRY[itype]
# end

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
