module Parser
using Configurations
using ArgCheck
using ..Utils
using ..Smearing

# Types and methods that are used by multiple parameter sets

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

@option struct QuollParams
    input::InputParams

    # Optional parameter sets
    output::Union{OutputParams, Nothing} = nothing
    basis_projection::Union{BasisProjectionParams, Nothing} = nothing

    # Optional parameter sets with existing defaults
    kpoint_grid::KPointGridParams = KPointGridParams()
    postprocessing::PostprocessParams = PostprocessParams()    
    error_metrics::ErrorMetricParams = ErrorMetricParams()
end

end