using Configurations
using ..Smearing

@option struct EigenvalueErrorParams <: AbstractQuollParams
    smearing_function::SmearingFunction = FermiDirac()
    temperature::Float64 = 1000.0
end

@option struct ElEntropyErrorParams <: AbstractQuollParams
    smearing_function::SmearingFunction = FermiDirac()
    temperature::Float64 = 1000.0
end

@option struct ErrorMetricParams <: AbstractQuollParams
    # TODO: could warn if reference or predicted data not found in supplied directories
    mae::Bool = false
    eigenvalue_error::Bool = false
    el_entropy_error::Bool = false
    eigenvalue_error_params::EigenvalueErrorParams = EigenvalueErrorParams()
    el_entropy_error_params::ElEntropyErrorParams = ElEntropyErrorParams()
end
