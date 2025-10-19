using AutoHashEquals
using Base

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
    n_basis_atom::Vector{Int}
    atom2species::Vector{ChemicalSpecies}
    basis2atom::Vector{Int}
    atom2basis::Vector{UnitRange{Int}}
end
# 1) Which basis functions belong to which atom? basis
# 2) which basis function INDICES belong to which atom? basis2atom
#    Do we need to know the indices in the case of every format?
#    It's something that's required if one wants to interconvert between certain formats.
#    More precisely, we need it when we convert it to a different format but we should
#    have freedom to change that indexing during conversions.

# If we load basis-indices and construct RealCSCSparsity, what do we get?
# We get a struct which describes the sparsity pattern a 3D array.
# Whereas basis2atom describes slices of that 3D array which belong
# to particular atoms.

# TODO: Alternatively (or additionally) could implement Base.getindex for Int and ChemicalSpecies
function basis_atom(basisset::BasisSetMetadata, iat::Int)
    return basisset.basis[basisset.atom2species[iat]]
end

function basis_species(basisset::BasisSetMetadata, species::ChemicalSpecies)
    return basisset.basis[species]
end

# Get a number of basis functions belonging to atoms with a lower index for every atom
function get_offsets(basisset::BasisSetMetadata)
    return [interval[begin] - 1 for interval in basisset.atom2basis]
end