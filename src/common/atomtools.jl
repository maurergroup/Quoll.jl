using AtomsIOPython

function recentre(atoms::AbstractSystem)
    # Currently the following function does not support recentering
    # for systems with open boundary conditions
    @argcheck all(periodicity(atoms))

    # Convert atom positions to fractional coordinates.
    # Strip units for matrix division because it doesn't always
    # work with units, resulting fractional coordinates do not have units
    # which is correct
    atom_pos = stack(position(atoms, :))
    lat = stack(cell_vectors(atoms))
    atom_frac_pos = ustrip(lat) \ ustrip(atom_pos)

    # Wrap fractional coordinates to [-0.5, 0.5) unit cell
    # and shift back to [0.0, 1.0)
    wrapped_atom_frac_pos = atom_frac_pos .- 1e-8
    wrapped_atom_frac_pos = wrapped_atom_frac_pos .- round.(wrapped_atom_frac_pos) .+ 1e-8
    center = SA[0.5, 0.5, 0.5]
    wrapped_atom_frac_pos .+= center

    # Remap fractional coordinates to atom positions
    wrapped_atom_pos = lat * wrapped_atom_frac_pos

    # Construct new atoms
    recentered_atoms = periodic_system(
        [species(atom) => wrapped_atom_pos[:, iat] for (iat, atom) in enumerate(atoms)],
        cell_vectors(atoms),
    )
    return recentered_atoms
end

# Given a dictionary z1z2... return an array ij...,
# leading to a dense data structure but with faster access (no hashing)
function speciesdict_to_atomarray(d, atom2species)
    firstkey = first(keys(d))
    firstkey isa ChemicalSpecies ? DIM = 1 : DIM = length(firstkey)

    concat = outer_concatenate(atom2species, Val(DIM))
    return get.(Ref(d), concat, missing)
end

outer_concatenate(vec, ::Val{1}) = vec
outer_concatenate(vec, ::Val{2}) = tuple.(vec, reshape(vec, 1, :))