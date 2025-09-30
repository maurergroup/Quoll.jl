module Parser
using Configurations
using ArgCheck
using ..Utils
using ..Smearing

abstract type AbstractQuollParams end

function Configurations.from_dict(
    ::Type{OT},
    ::OptionField{:smearing_function},
    ::Type{Type{<:SmearingFunction}},
    s
) where OT<:AbstractQuollParams
    @argcheck isa(s, String)
    
    symbol = Symbol(Utils.normalize_comparison(s))
    @argcheck hasmethod(implemented_smearing, Tuple{Val{symbol}}) "Smearing function $s is unavailable or doesn't exist"

    return implemented_smearing(Val(symbol))
end

include("Parser/input_output.jl")
include("Parser/postprocessing.jl")
include("Parser/basis_projection.jl")
include("Parser/error_metrics.jl")
include("Parser/kpoint_grid.jl")

"""
    read_inputfile(filepath::AbstractString)

Read Quoll input file in TOML format.
"""
function read_inputfile(filepath::AbstractString) end

# @option struct QuollParams
#     input::InputParams
#     output::Union{OutputParams, Nothing} = nothing
# end

end