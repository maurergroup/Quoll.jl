module Basis
using Dictionaries
using AtomsBase
using AutoHashEquals

export BasisMetadata, BasisSetMetadata

@auto_hash_equals struct BasisMetadata{E}
    z::ChemicalSpecies
    n::Int
    l::Int
    m::Int
    extras::E
end

struct BasisSetMetadata{E}
    basis::Dictionary{ChemicalSpecies,Vector{BasisMetadata{E}}}
    n_basis_atom::Vector{Int}
    atom2species::Vector{ChemicalSpecies}
    basis2atom::Vector{Int}
    atom2basis::Vector{UnitRange{Int}}
end

end