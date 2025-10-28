using Base
using AutoHashEquals

@auto_hash_equals struct BasisMetadata{E}
    z::ChemicalSpecies
    n::Int
    l::Int
    m::Int
    extras::E
end

function Base.show(io::IO, basisf::BasisMetadata)
    print(io, "$(basisf.z)$(basisf.n)$(basisf.l)$(basisf.m)")
end

# Would BasisSetMetadata have to be modified in the case of SOC?
# I think a reasonable (although maybe not the most memory-efficient) way would be to assume
# BasisSetMetadata.basis would store 2x basis functions without realising some of them are
# spin up and spin down. However, can we tell which basis function would have which spin?
# That would depend on electronic structure method (in which case SpinsMetadata constructor
# could be additionally dispatched on Type{AbstractOperator}).
struct BasisSetMetadata{E}
    basis::Dictionary{ChemicalSpecies, Vector{BasisMetadata{E}}}
    atom2species::Vector{ChemicalSpecies}
end
# 1) Which basis functions belong to which atom? basis
# 2) which basis function INDICES belong to which atom? basis2atom
#    Do we need to know the indices in the case of every format?
#    It's something that's required if one wants to interconvert between certain formats.
#    More precisely, we need it when we convert it to a different format but we should
#    have freedom to change that indexing during conversions.

# If we load basis-indices and construct RealCSCSparsity, what do we get?
# We get a struct which describes the sparsity pattern of a 3D array.
# Whereas basis2atom describes regions of that 3D array which belong
# to particular atoms.


get_atom2nbasis(basisset::BasisSetMetadata) = get_atom2nbasis(basisset.basis, basisset.atom2species)
get_atom2nbasis(basis, atom2species) = [length(basis[z]) for z in atom2species]

get_atom2basis(basisset::BasisSetMetadata) = get_atom2basis(basisset.basis, basisset.atom2species)

function get_atom2basis(basis, atom2species)
    atom2nbasis = get_atom2nbasis(basis, atom2species)
    return get_atom2basis(atom2nbasis)
end

function get_atom2basis(atom2nbasis)
    n_atoms = length(atom2nbasis)
    atom2basis_start = cumsum(atom2nbasis) .- atom2nbasis .+ 1
    return [
        atom2basis_start[iat]:atom2basis_start[iat] + atom2nbasis[iat] - 1
        for iat in 1:n_atoms
    ]
end

get_basis2atom(basisset::BasisSetMetadata) = get_basis2atom(basisset.basis, basisset.atom2species)

function get_basis2atom(basis, atom2species)
    atom2nbasis = get_atom2nbasis(basis, atom2species)
    return get_basis2atom(atom2nbasis)
end

get_basis2atom(atom2nbasis) = reduce(vcat, [fill(iat, nbasis) for (iat, nbasis) in enumerate(atom2nbasis)])

# Get a number of basis functions for a given chemical species
function get_species2nbasis(basisset::BasisSetMetadata)
    return dictionary([
        z => length(basis_species(basisset, z))
        for z in keys(basisset.basis)
    ])
end

# TODO: Alternatively (or additionally) could implement Base.getindex for Int and ChemicalSpecies
function basis_atom(basisset::BasisSetMetadata, iat::Int)
    return basisset.basis[basisset.atom2species[iat]]
end

function basis_species(basisset::BasisSetMetadata, species::ChemicalSpecies)
    return basisset.basis[species]
end

# Get a number of basis functions belonging to atoms with a lower index for every atom
function get_offsets(basisset::BasisSetMetadata)
    atom2basis = get_atom2basis(basisset)
    return get_offsets(atom2basis)
end

get_offsets(atom2basis) = [interval[begin] - 1 for interval in atom2basis]
    