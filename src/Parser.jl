module Parser
using Configurations
using ArgCheck
using ..Utils
using ..Smearing

abstract type AbstractQuollParams end

function Configurations.from_dict(
    ::Type{OT},
    ::OptionField{:temperature},
    ::Type{Float64},
    temperature
) where OT<:AbstractQuollParams
    if isa(temperature, String)
        temperature = parse(Float64, temperature)
    elseif !(isa(temperature, Float64))
        temperature = convert(Float64, temperature)
    end
    @argcheck temperature ≥ 0.0
    return temperature
end

function Configurations.from_dict(
    ::Type{OT},
    ::OptionField{:smearing_function},
    ::Type{SmearingFunction},
    smearing_function
) where OT<:AbstractQuollParams
    @argcheck isa(smearing_function, String)

    itype = findfirst(
        x -> x == Utils.normalize_comparison(smearing_function),
        lowercase.(string.(nameof.(SMEAR_REGISTRY)))
    )
    @argcheck itype !== nothing "Smearing function $smearing_function not found"

    return SMEAR_REGISTRY[itype]()
end

include("Parser/input_output.jl")
include("Parser/postprocessing.jl")
include("Parser/basis_projection.jl")
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