module Parser

using Configurations
using StaticArrays
using DelimitedFiles
using Unitful
using UnitfulAtomic
using AtomsBase
using ArgCheck
using MPI

using ..Quoll
using ..Quoll:
    AbstractOperator,
    OperatorKind,
    KGrid,
    allows_symmetry,
    get_avail_operatorkinds,
    get_operatorkinds,
    get_readformat,
    get_writeformat
import ..Quoll: find_operatorkinds, get_kgrid

using ..Projections
using ..Projections: AbstractBasisProjection, get_basis_projection
import ..Projections: perform_core_projection

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
include("parser/symmetry.jl")

@option struct QuollParams
    input::InputParams

    # Optional parameter sets
    output::Maybe{OutputParams} = nothing
    basis_projection::Maybe{BasisProjectionParams} = nothing

    # Optional parameter sets with existing defaults
    kpoint_grid::KPointGridParams = KPointGridParams()
    postprocessing::PostprocessParams = PostprocessParams()
    error_metrics::ErrorMetricParams = ErrorMetricParams()
    symmetry::SymmetryParams = SymmetryParams()

    function QuollParams(
        input,
        output,
        basis_projection,
        kpoint_grid,
        postprocessing,
        error_metrics,
        symmetry,
    )
        # Check if input operators are appropriate for the task
        validate_operatorkinds(
            input.operators,
            output,
            basis_projection,
            postprocessing,
            error_metrics,
        )
        # Check if supplied symmetry options are appropriate
        validate_symmetry(input, basis_projection, symmetry)

        # Look if there are any unresolvable clashes
        search_clashes(basis_projection, error_metrics)

        new(input, output, basis_projection, kpoint_grid, postprocessing, error_metrics, symmetry)
    end
end

# Methods that either require subparameter types or QuollParams itself
include("parser/methods.jl")

end