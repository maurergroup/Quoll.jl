module Basis
using Dictionaries
using AtomsBase
using AutoHashEquals
using Dictionaries
using ArgCheck

using ..OperatorIO

export BasisMetadata, BasisSetMetadata

@auto_hash_equals struct BasisMetadata{E}
    z::ChemicalSpecies
    n::Int
    l::Int
    m::Int
    extras::E
end

struct BasisSetMetadata{E}
    basis::Dictionary{ChemicalSpecies, Vector{BasisMetadata{E}}}
    n_basis_atom::Vector{Int}
    atom2species::Vector{ChemicalSpecies}
    basis2atom::Vector{Int}
    atom2basis::Vector{UnitRange{Int}}
end

function BasisSetMetadata(dir::AbstractString, atoms::AbstractSystem, metadata::FHIaimsOperatorMetadata)
    p = joinpath(dir, "basis-indices.out")
    @argcheck ispath(p)

    n_atoms = length(atoms)
    n_basis_atom = zeros(Int, n_atoms)
    atom2species = species(atoms, :)
    basis2atom = Vector{Int}(undef, metadata.n_basis)

    basis = dictionary(
        species => BasisMetadata{Dict{String, String}}[]
        for species in unique(atom2species)
    )

    f = open(p, "r")

    # Reads both the empty lines and the header
    while isempty(strip(readline(f))) end

    species2firstatom = dictionary(
        species => findfirst(x -> x == species, atom2species)
        for species in unique(atom2species)
    )
    for ib in 1:metadata.n_basis
        _, type, iat, n, l, m = convert.(String, split(readline(f)))
        iat, n, l, m = parse.(Int, (iat, n, l, m))

        if iat == species2firstatom[atom2species[iat]]
            push!(
                basis[atom2species[iat]],
                BasisMetadata(atom2species[iat], n, l, m, Dict{String, String}("type" => type))
            )
        end

        basis2atom[ib] = iat
        n_basis_atom[iat] += 1
    end
    close(f)

    atom2basis_start = cumsum(n_basis_atom) .- n_basis_atom .+ 1
    atom2basis = [
        atom2basis_start[iat]:atom2basis_start[iat] + n_basis_atom[iat] - 1
        for iat in 1:n_atoms
    ]

    return BasisSetMetadata(basis, n_basis_atom, atom2species, basis2atom, atom2basis)
end

end