using Configurations
using StaticArrays
using ..Smearing

"""
"""
struct KPathParams
    kpoint_begin::SVector{3, Float64}
    kpoint_end::SVector{3, Float64}
    n_points::Int
    symbol_begin::Symbol
    symbol_end::Symbol
end

@option struct FermiLevelParams <: AbstractQuollParams
    smearing_function::SmearingFunction = Gaussian()
end

@option struct DOSParams <: AbstractQuollParams
    smearing_function::SmearingFunction = Gaussian()
    temperature::Float64 = 1000.0
    n_points::Int = 1000
    energy_begin::Float64 = -5.0
    energy_end::Float64 = 5.0

    function DOSParams(smearing_function, temperature, n_points, energy_begin, energy_end)
        @argcheck energy_begin < energy_end
        new(smearing_function, temperature, n_points, energy_begin, energy_end)
    end
end

@option struct BandStructureParams <: AbstractQuollParams
    # TODO: Should convert bands to a KPath from Brillouin
    # so that I can use utilities like interpolate and plotting
    # once I have atoms information
    bands::Union{Vector{KPathParams}, Nothing} = nothing
end

@option struct PostprocessParams <: AbstractQuollParams
    fermi_level::Bool = true
    dos::Bool = true
    band_structure::Bool = true
    fermi_level_params::FermiLevelParams = FermiLevelParams()
    dos_params::DOSParams = DOSParams()
    band_structure_params::BandStructureParams = BandStructureParams()
end

function parse_kpathparams(band::AbstractString)
    matches = getfield.(eachmatch(r"\S+", band), :match)
    @argcheck length(matches) == 9

    kpoint_begin = tryparse.(Float64, matches[1:3])
    kpoint_end = tryparse.(Float64, matches[4:6])
    n_points = tryparse(Int, matches[7])

    @argcheck all(x -> x !== nothing, kpoint_begin) && all(x -> x !== nothing, kpoint_end) && n_points !== nothing "K-point values in the band structure are not parsable"

    symbol_begin = Symbol(matches[8])
    symbol_end = Symbol(matches[9])

    return KPathParams(SVector{3, Float64}(kpoint_begin), SVector{3, Float64}(kpoint_end), n_points, symbol_begin, symbol_end)
end

function Configurations.from_dict(
    ::Type{BandStructureParams},
    ::OptionField{:bands},
    ::Type{Vector{KPathParams}},
    bands,
)
    @argcheck isa(bands, Vector{String})

    # Vector{String} -> Vector{KPathParams}
    return parse_kpathparams.(bands)
end