module Basis
using Dictionaries
using AtomsBase
export BasisMetadata, BasisSetMetadata

struct BasisMetadata
    z::ChemicalSpecies
    n::Int
    l::Int
    m::Int
    extra::NamedTuple
end

struct BasisSetMetadata
    basis::Dictionary{ChemicalSpecies,Vector{BasisMetadata}}
    n_basis_atom::Vector{Int}
    atom2species::Vector{ChemicalSpecies}
    basis2atom::Vector{Int}
    atom2basis::Vector{UnitRange{Int}}
end

end