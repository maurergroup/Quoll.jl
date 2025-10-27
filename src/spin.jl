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

# Constructors for SpinsMetadata from OperatorKind, BasisSetMetadata, and AbstractOperator.
# AbstractOperator is required for the SOC case. This is because for now in the case of SOC
# BasisSetMetadata.basis stores 2x the number of basis functions, and we don't know which of them
# are spin up and spin down. This type of information might depend on the format

function SpinsMetadata(kind::OperatorKind, basisset::BasisSetMetadata, format::Type{<:AbstractOperator})
    return SpinsMetadata(kind, basisset.basis, format)
end

function SpinsMetadata(kind::Hamiltonian, basis::Dictionary{ChemicalSpecies}, format::Type{<:AbstractOperator})
    return SpinsMetadata(kind, Val(kind.spin), basis, format)
end

function SpinsMetadata(::Hamiltonian, ::Val{:up}, basis::Dictionary{ChemicalSpecies}, ::Type{<:AbstractOperator})
    return SpinsMetadata(
        Dictionary(
            keys(basis),
            [fill(⬆, length(atom_basis)) for atom_basis in basis]
        ),
        false,
    )
end

function SpinsMetadata(::Hamiltonian, ::Val{:down}, basis::Dictionary{ChemicalSpecies}, ::Type{<:AbstractOperator})
    return SpinsMetadata(
        Dictionary(
            keys(basis),
            [fill(⬇, length(atom_basis)) for atom_basis in basis]
        ),
        false,
    )
end

SpinsMetadata(::Hamiltonian, ::Val{:none}, ::Dictionary{ChemicalSpecies}, ::Type{<:AbstractOperator}) = nothing

SpinsMetadata(::Overlap, ::Dictionary{ChemicalSpecies}, ::Type{<:AbstractOperator}) = nothing