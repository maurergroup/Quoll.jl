module Parser

using Configurations
using StaticArrays
using Unitful
using UnitfulAtomic
using ArgCheck

using ..Projections
using ..Postprocessing

# Common utilities used in multiple instances
include("parser/common.jl")

# Subparameter types and methods to parse them
include("parser/input_output.jl")
include("parser/postprocessing.jl")
include("parser/basis_projection.jl")
include("parser/error_metrics.jl")
include("parser/kpoint_grid.jl")

@option struct QuollParams
    input::InputParams

    # Optional parameter sets
    output::Union{OutputParams, Nothing} = nothing
    basis_projection::Union{BasisProjectionParams, Nothing} = nothing

    # Optional parameter sets with existing defaults
    kpoint_grid::KPointGridParams = KPointGridParams()
    postprocessing::PostprocessParams = PostprocessParams()    
    error_metrics::ErrorMetricParams = ErrorMetricParams()

    function QuollParams(input, output, basis_projection, kpoint_grid, postprocessing, error_metrics)
        # Check if input operators are appropriate for the task
        validate_operatorkinds(input.operators, output, basis_projection, postprocessing, error_metrics)
        new(input, output, basis_projection, kpoint_grid, postprocessing, error_metrics)
    end
end

# Methods that either require all subparameter types or QuollParams itself
include("parser/validation.jl")

end