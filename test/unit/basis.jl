# TODO: write unit tests

# function basis_atom(basisset::BasisSetMetadata, iat::Int)
#     return basisset.basis[basisset.atom2species[iat]]
# end

# function basis_species(basisset::BasisSetMetadata, species::ChemicalSpecies)
#     return basisset.basis[species]
# end

# # Get a number of basis functions belonging to atoms with a lower index for every atom
# function get_offsets(basisset::BasisSetMetadata)
#     return [interval[begin] - 1 for interval in basisset.atom2basis]
# end
    
# # Get a number of basis functions for a given chemical species
# function get_species2nb(basisset::BasisSetMetadata)
#     return dictionary([
#         z => length(basis_species(basisset, z))
#         for z in keys(basisset.basis)
#     ])
# end