module Parser

using Configurations
using StaticArrays
using DelimitedFiles
using Unitful
using UnitfulAtomic
using AtomsBase
using ArgCheck

using ..Quoll
using ..Quoll:
    AbstractOperator,
    OperatorKind,
    get_avail_operatorkinds,
    get_operatorkinds,
    get_readformat,
    get_writeformat
import ..Quoll: find_operatorkinds
using ..Projections
using ..Projections: AbstractBasisProjection, get_basis_projection
using ..Postprocessing
using ..Postprocessing: SmearingFunction, get_smearing

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
    output::Union{OutputParams,Nothing} = nothing
    basis_projection::Union{BasisProjectionParams,Nothing} = nothing

    # Optional parameter sets with existing defaults
    kpoint_grid::KPointGridParams = KPointGridParams()
    postprocessing::PostprocessParams = PostprocessParams()
    error_metrics::ErrorMetricParams = ErrorMetricParams()

    function QuollParams(
        input,
        output,
        basis_projection,
        kpoint_grid,
        postprocessing,
        error_metrics,
    )
        # Check if input operators are appropriate for the task
        validate_operatorkinds(
            input.operators,
            output,
            basis_projection,
            postprocessing,
            error_metrics,
        )
        search_clashes(basis_projection, error_metrics)
        new(input, output, basis_projection, kpoint_grid, postprocessing, error_metrics)
    end
end

# Methods that either require all subparameter types or QuollParams itself
include("parser/methods.jl")

end