using Configurations
using ..Smearing

@option struct EigenvalueErrorParams <: AbstractQuollParams
    smearing_function::SmearingFunction = FermiDirac()
    temperature::Float64 = 1000.0

    function EigenvalueErrorParams(smearing_function, temperature)
        @argcheck temperature > 0.0
        new(smearing_function, temperature)
    end
end

@option struct ElEntropyErrorParams <: AbstractQuollParams
    smearing_function::SmearingFunction = FermiDirac()
    temperature::Float64 = 1000.0

    function ElEntropyErrorParams(smearing_function, temperature)
        @argcheck temperature > 0.0
        new(smearing_function, temperature)
    end
end

@option struct ErrorMetricParams <: AbstractQuollParams
    # TODO: could warn if reference or predicted data not found in supplied directories
    # TODO: Also I could allow for subblock MAE here
    mae::Bool = false
    eigenvalue_error::Bool = false
    el_entropy_error::Bool = false
    eigenvalue_error_params::EigenvalueErrorParams = EigenvalueErrorParams()
    el_entropy_error_params::ElEntropyErrorParams = ElEntropyErrorParams()
end
