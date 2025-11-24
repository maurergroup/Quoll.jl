function get_spglib_cell(atoms::FlexibleSystem)
    real_lattice = get_real_lattice(atoms)
    real_frac_positions = get_frac_positions(atoms)
    atomic_numbers = atomic_number(atoms, :)

    return Cell(ustrip.(real_lattice), real_frac_positions, atomic_numbers)
end

function get_symmetry_rotations(spglib_cell; crystal_symmetry = false, symprec = 1e-5)
    if crystal_symmetry
        rotations, _ = get_symmetry(spglib_cell, symprec)
    else
        rotations = [SMatrix{3, 3}(I(3))]
    end
    return rotations
end
