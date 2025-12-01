@enum Spin begin
    ⬆ = 1
    ⬇ = -1
end

struct SpinsMetadata{S<:SpeciesAnyDict}
    spins::S
    soc::Bool
end

function SpinsMetadata(spins::SpinsMetadata{<:SpeciesDictionary}, masks::SpeciesDictionary{BitVector})
    subspins_dict = getindex.(spins.spins, masks)
    return SpinsMetadata(subspins_dict, spins.soc)
end

# TODO: what about atom2species?
function spins_species(spins::SpinsMetadata, z::ChemicalSpecies)
    return spins.spins[z]
end
