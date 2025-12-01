function recentre(atoms::AbstractSystem)
    # Currently the following function does not support recentering
    # for systems with open boundary conditions
    @argcheck all(periodicity(atoms))

    # Convert atom positions to fractional coordinates.
    atom_frac_pos = get_frac_positions(atoms)

    # Wrap fractional coordinates to [-0.5, 0.5) unit cell
    # and shift back to [0.0, 1.0)
    wrapped_atom_frac_pos = atom_frac_pos .- 1e-8
    wrapped_atom_frac_pos = wrapped_atom_frac_pos .- round.(wrapped_atom_frac_pos) .+ 1e-8
    center = SA[0.5, 0.5, 0.5]
    wrapped_atom_frac_pos .+= center

    # Remap fractional coordinates to atom positions
    wrapped_atom_pos = get_real_lattice(atoms) * wrapped_atom_frac_pos

    # Construct new atoms
    recentered_atoms = periodic_system(
        [Atom(atom, position = wrapped_atom_pos[:, iat]) for (iat, atom) in enumerate(atoms)],
        cell_vectors(atoms),
    )
    return recentered_atoms
end

function get_real_lattice(atoms::AbstractSystem)
    return stack(cell_vectors(atoms))
end

function get_recip_lattice(atoms::AbstractSystem)
    real_lattice = get_real_lattice(atoms)
    return get_recip_lattice(real_lattice)
end

# A^T B = 2πI
# where both matrices contain basis vectors as columns
function get_recip_lattice(real_lattice)
    return transpose(ustrip.(real_lattice)) \ (2π * I(3))
end

# Strip units for matrix division because it doesn't always work with units.
# Resulting fractional coordinates do not have units which is correct
function get_frac_positions(atoms::AbstractSystem)
    real_lattice = get_real_lattice(atoms)
    real_positions = ustrip.(stack(position(atoms, :)))
    return ustrip.(real_lattice) \ ustrip.(real_positions)
end

function get_species2atom(atoms::AbstractSystem)
    all_species = species(atoms, :)
    unique_species = Indices(unique(all_species))
    return map(unique_species) do z
        findall(isequal(z), all_species)
    end
end