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
        Sref = Overlap(:ref) ∈ input.operators
        Spred = Overlap(:pred) ∈ input.operators
        Href = Hamiltonian(:ref) ∈ input.operators
        Hpred = Hamiltonian(:pred) ∈ input.operators

        output !== nothing && (@argcheck Sref || Spred || Href || Hpred)
        basis_projection !== nothing && (@argcheck Sref)

        postprocessing.dos && (@argcheck (Sref && Href) || (Spred && Hpred))
        postprocessing.fermi_level && (@argcheck (Sref && Href) || (Spred && Hpred))
        postprocessing.band_structure && (@argcheck (Sref && Href) || (Spred && Hpred))

        error_metrics.mae && (@argcheck Href && Hpred)
        error_metrics.eigenvalue_error && (@argcheck (Sref && Href) || (Spred && Hpred))
        error_metrics.el_entropy_error && (@argcheck (Sref && Href) || (Spred && Hpred))

        new(input, output, basis_projection, kpoint_grid, postprocessing, error_metrics)
    end
end

# TODO: remove
# output               Sref || Spred || Href || Hpred
# basis_projection     Sref
# kpoint_grid          Sref || Spred || Href || Hpred
# postprocessing
#   dos                (Sref && Href) || (Spred && Hpred)
#   fermi_level        (Sref && Href) || (Spred && Hpred)
#   band_structure     (Sref && Href) || (Spred && Hpred)
# error_metrics
#   mae                Href && Hpred
#   eigenvalue_error   (Sref && Href) || (Spred && Hpred)
#   el_entropy_error   (Sref && Href) || (Spred && Hpred)

end