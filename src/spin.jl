using Base

@enum Spin begin
    ⬆ = 1
    ⬇ = -1
end

struct SpinsMetadata{S<:SpeciesAnyDict}
    spins::S
    soc::Bool
end

# TODO: what about atom2species?
function spins_species(spins::SpinsMetadata, z::ChemicalSpecies)
    return spins.spins[z]
end
