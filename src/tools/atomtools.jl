"""
    recentre(atoms) -> AbstractSystem

Wrap atom positions into the unit cell centred at [0.5, 0.5, 0.5] in fractional
coordinates. Only supports fully periodic systems. Returns a new system with recentred
positions.
"""
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
        [Atom(atom; position=wrapped_atom_pos[:, iat]) for (iat, atom) in enumerate(atoms)],
        cell_vectors(atoms),
    )
    return recentered_atoms
end

"""
    get_real_lattice(atoms) -> Matrix

Return the real-space lattice as a matrix with basis vectors as columns.
"""
function get_real_lattice(atoms::AbstractSystem)
    return stack(cell_vectors(atoms))
end

"""
    get_recip_lattice(atoms_or_real_lattice) -> Matrix

Return the reciprocal-space lattice satisfying `AᵀB = 2πI`, with basis vectors as columns.
"""
function get_recip_lattice(atoms::AbstractSystem)
    real_lattice = get_real_lattice(atoms)
    return get_recip_lattice(real_lattice)
end

# A^T B = 2πI
# where both matrices contain basis vectors as columns
function get_recip_lattice(real_lattice)
    return transpose(ustrip.(real_lattice)) \ (2π * I(3))
end

"""
    get_frac_positions(atoms) -> Matrix

Return atom positions in fractional coordinates (unitless). Each column is one atom.
"""
function get_frac_positions(atoms::AbstractSystem)
    real_lattice = get_real_lattice(atoms)
    real_positions = ustrip.(stack(position(atoms, :)))
    return ustrip.(real_lattice) \ ustrip.(real_positions)
end

"""
    get_species2atom(atoms) -> Dictionary

Return a dictionary mapping each unique `ChemicalSpecies` to the vector of atom indices of
that species.
"""
function get_species2atom(atoms::AbstractSystem)
    all_species = species(atoms, :)
    unique_species = Indices(unique(all_species))
    return map(unique_species) do z
        findall(isequal(z), all_species)
    end
end