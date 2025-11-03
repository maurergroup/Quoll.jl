using Base

@enum Spin begin
    ⬆ = 1
    ⬇ = -1
end

struct SpinsMetadata{S<:SpeciesDict{Spin, 1}}
    spins::S
    soc::Bool
end

# TODO: what about atoms2species?
function spins_species(spins::SpinsMetadata, z::ChemicalSpecies)
    return spins.spins[z]
end
