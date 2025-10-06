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
    @argcheck hasmethod(smearing, Tuple{Val{symbol}}) "Smearing function $s is unavailable or doesn't exist"

    return smearing(Val(symbol))
end

include("Parser/input_output.jl")
include("Parser/postprocessing.jl")
include("Parser/basis_projection.jl")
include("Parser/error_metrics.jl")
include("Parser/kpoint_grid.jl")

# Check if the operators are sufficient for requested job
function validate_operatorkinds(operators, output, basis_projection, postprocessing, error_metrics)
    Sref = Overlap(:ref) ∈ operators
    Spred = Overlap(:pred) ∈ operators
    Href = Hamiltonian(:ref) ∈ operators
    Hpred = Hamiltonian(:pred) ∈ operators

    output !== nothing && (@argcheck Sref || Spred || Href || Hpred)
    basis_projection !== nothing && (@argcheck Sref)

    postprocessing.dos && (@argcheck (Sref && Href) || (Spred && Hpred))
    postprocessing.fermi_level && (@argcheck (Sref && Href) || (Spred && Hpred))
    postprocessing.band_structure && (@argcheck (Sref && Href) || (Spred && Hpred))

    error_metrics.mae && (@argcheck Href && Hpred)
    error_metrics.eigenvalue_error && (@argcheck (Sref && Href) || (Spred && Hpred))
    error_metrics.el_entropy_error && (@argcheck (Sref && Href) || (Spred && Hpred))

    return true
end

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

end