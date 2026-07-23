"""
    recentre(atoms) -> AbstractSystem

Wrap atom positions to the (-0.5, 0.5] fractional coordinate range. Only supports fully
periodic systems. Returns a new system with recentred positions.
"""
function recentre(atoms::AbstractSystem)
    # Currently the following function does not support recentering
    # for systems with open boundary conditions
    @argcheck all(periodicity(atoms))

    # Convert atom positions to fractional coordinates.
    atom_frac_pos = get_frac_positions(atoms)

    # Wrap fractional coordinates to [-0.5, 0.5) unit cell
    wrapped_atom_frac_pos = atom_frac_pos .- 1e-8
    wrapped_atom_frac_pos = wrapped_atom_frac_pos .- round.(wrapped_atom_frac_pos) .+ 1e-8

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

"""
    get_images(atoms) -> Matrix

Return periodic images to which each atom belongs. Each column is one atom.
"""
function get_images(atoms)
    frac_positions = Quoll.get_frac_positions(atoms)
    images = floor.(Int, frac_positions)
    return [SVector{3}(img) for img in eachcol(images)]
end

"""
    get_centroid(coords) -> Matrix

Return the centroid (mean position) of a `3 × N` coordinate matrix as a `3 × 1` matrix.
"""
function get_centroid(coords::AbstractMatrix)
    norm = 1 / size(coords, 2)
    return norm * sum(coords; dims=2)
end

"""
    get_translation_invariant(coords) -> Matrix

Return a translation-invariant version of a `3 × N` coordinate matrix by subtracting its
[`get_centroid`](@ref). Two configurations that differ only by a rigid translation compare
equal after this transformation.
"""
function get_translation_invariant(coords::AbstractMatrix)
    centroid = get_centroid(coords)
    return coords .- centroid
end

"""
    same_atoms(atoms1, atoms2) -> Bool

Return `true` if two atomic systems are equivalent for the purpose of shared atom indexing:
same chemical species (in the same order), the same cell, and the same atom positions up to a
rigid translation (absolute positions may differ, only the translation-invariant positions are
compared). Positions and cell vectors are assumed to share a common length unit.
"""
function same_atoms(atoms1::AbstractSystem, atoms2::AbstractSystem)
    # Same species in the same order (also catches a differing number of atoms).
    species(atoms1, :) == species(atoms2, :) || return false

    # Same cell.
    cell1 = ustrip.(get_real_lattice(atoms1))
    cell2 = ustrip.(get_real_lattice(atoms2))
    isapprox(cell1, cell2) || return false

    # Same atom positions up to a rigid translation.
    pos1 = get_translation_invariant(ustrip.(stack(position(atoms1, :))))
    pos2 = get_translation_invariant(ustrip.(stack(position(atoms2, :))))
    return isapprox(pos1, pos2)
end
