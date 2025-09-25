using Configurations

@option struct FermiLevelParams
    smearing_function::String = "fermi_dirac"
end

@option struct DOSParams
    smearing_function::String = "fermi_dirac"
    smearing::Float64 = 0.1
    n_points::Int = 1000
    start_energy::Float64
    end_energy::Float64
end

@option struct BandStructureParams
    # TODO: Should be something else compatible with Brillouin
    bands::Vector{String}
end

@option struct PostprocessParams
    fermi_level::Bool = true
    fermi_level_params::FermiLevelParams
    dos::Bool = true
    dos_params::DOSParams
    band_structure::Bool = true
    band_structure_params::BandStructureParams
end
