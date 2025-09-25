module Parser
using Configurations

include("Parser/input_output.jl")
include("Parser/postprocessing.jl")
include("Parser/core_projection.jl")
include("Parser/error_metrics.jl")

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