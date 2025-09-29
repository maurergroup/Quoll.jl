using Configurations

@option struct FermiLevelParams <: AbstractQuollParams
    smearing_function::String = "fermi_dirac"
end

@option struct DOSParams <: AbstractQuollParams
    smearing_function::String = "fermi_dirac"
    smearing::Float64 = 0.1
    n_points::Int = 1000
    start_energy::Float64
    end_energy::Float64
end

@option struct BandStructureParams <: AbstractQuollParams
    # TODO: Should be something else compatible with Brillouin
    bands::Vector{String}
end

@option struct PostprocessParams <: AbstractQuollParams
    fermi_level::Bool = true
    fermi_level_params::FermiLevelParams
    dos::Bool = true
    dos_params::DOSParams
    band_structure::Bool = true
    band_structure_params::BandStructureParams
end
