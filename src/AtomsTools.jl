module AtomsTools
using AtomsIOPython
using AtomsBase
using StaticArrays
using Unitful
using ArgCheck
using ..OperatorIO
export recentre, load_atoms

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

function load_atoms(dir::AbstractString, ::Type{FHIaimsOperator})
    path = joinpath(dir, "geometry.in")
    atoms = load_system(path)

    # Recentre atom positions to be consistent with relative positions
    # that are used during real-space matrix integration in FHI-aims
    all(periodicity(atoms)) && (atoms = recentre(atoms))
    
    return atoms
end

end