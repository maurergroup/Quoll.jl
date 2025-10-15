using Base

@enum Spin begin
    ⬆ = 1
    ⬇ = -1
end

struct SpinsMetadata
    spins::Dictionary{ChemicalSpecies, Vector{Spin}}
    soc::Bool
end
