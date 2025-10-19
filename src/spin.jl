using Base

@enum Spin begin
    ⬆ = 1
    ⬇ = -1
end

struct SpinsMetadata
    spins::Dictionary{ChemicalSpecies, Vector{Spin}}
    soc::Bool
end

# TODO: what about atoms2species?
function spins_species(spins::SpinsMetadata, z::ChemicalSpecies)
    return spins.spins[z]
end
