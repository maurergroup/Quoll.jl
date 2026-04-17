"""
    KGrid{T,W}

Container for an irreducible k-point grid with integration weights.
    
Fields:
- `kpoints::T` — collection of irreducible k-point coordinate vectors (fractional).
- `weights::W` — integration weights summing to 1.
- `time_reversal::Bool` — whether time-reversal symmetry was used to reduce the grid.
- `crystal_symmetry::Bool` — whether crystal point-group symmetry was used to reduce the grid.
"""
struct KGrid{T,W}
    kpoints::T
    weights::W
    time_reversal::Bool
    crystal_symmetry::Bool
end

get_nkpoints(kgrid::KGrid) = length(kgrid.kpoints)
get_kpoints(kgrid::KGrid, ik) = kgrid.kpoints[ik]
get_weights(kgrid::KGrid, ik) = kgrid.weights[ik]

"""
    construct_kgrid(atoms; density, mesh, shift, time_reversal, crystal_symmetry, symprec)

Build an irreducible k-point grid for the given atomic system using Spglib.

Returns a [`KGrid`](@ref) with irreducible k-points and integration weights. The grid
is reduced using the requested symmetries and validated by checking that all reducible
k-points can be regenerated. Falls back to no symmetries if validation fails.

The function works only for periodic, or fully non-periodic systems (in which case a gamma
point is returned).

The implementation is inspired by [DFTK.jl](https://github.com/JuliaMolSim/DFTK.jl).

# Keyword arguments
- `density=nothing` — reciprocal-space density (points per Å⁻¹ per axis). Ignored if `mesh`
  is given.
- `mesh=nothing` — explicit Monkhorst–Pack grid dimensions `[n₁, n₂, n₃]`. Takes priority
  over `density`.
- `shift=falses(3)` — half-grid shift per axis (`true` shifts by half a grid spacing).
- `time_reversal=false` — exploit time-reversal symmetry (k ↔ −k) to further reduce the grid.
- `crystal_symmetry=false` — exploit crystal point-group symmetry to reduce the grid.
- `symprec=1e-5` — symmetry detection tolerance (Å) passed to Spglib.
"""
function construct_kgrid(
    atoms::AbstractSystem;
    density=nothing,
    mesh=nothing,
    shift=falses(3),
    time_reversal=false,
    crystal_symmetry=false,
    symprec=1e-5,
)
    # Check the periodicity of the system
    if !any(periodicity(atoms))
        return KGrid([SVector{3}([0.0, 0.0, 0.0])], [1.0], false, false)
    elseif !all(periodicity(atoms))
        throw(error("`construct_kgrid` does not support systems with mixed boundary conditions"))
    end

    # Get symmetry rotations
    spglib_cell = get_spglib_cell(atoms)
    rotations = get_symmetry_rotations(
        spglib_cell; crystal_symmetry=crystal_symmetry, symprec=symprec
    )

    # Get k-point grid
    mesh = get_recip_mesh(atoms, density, mesh)
    kgrid_spglib = get_stabilized_reciprocal_mesh(
        rotations, mesh; is_shift=shift, is_time_reversal=time_reversal
    )

    # Add the shift and wrap irreducible k-points to [-0.5, 0.5)
    kirreds = map(Spglib.eachpoint(kgrid_spglib)) do kcoord
        normalize_kpoint_coordinate(SVector{3}(kcoord))
    end

    # Get groups of reducible k-points
    ir_mapping = kgrid_spglib.ir_mapping_table
    k_all_reducible = [findall(isequal(elem), ir_mapping) for elem in unique(ir_mapping)]

    # Get k-point weights
    n_equivalent_k = length.(k_all_reducible)
    weights = n_equivalent_k / sum(n_equivalent_k)

    # Test if all reducible k-points can be reproduced with symmetry rotations
    grid_address = kgrid_spglib.grid_address
    symmetries_valid = check_kpoint_reduction(
        rotations, time_reversal, shift, mesh, grid_address, k_all_reducible, kirreds
    )

    if !symmetries_valid
        @warn "Reducible k-points could not be generated from the irreducible kpoints. " *
            "This points to a bug in spglib. Defaulting to no symmetries."
        return construct_kgrid(
            atoms; density=density, mesh=mesh, shift=shift,
            time_reversal=false, crystal_symmetry=false, symprec=symprec,
        )
    else
        return KGrid(
            collect(eachpoint(kgrid_spglib)), weights, time_reversal, crystal_symmetry
        )
    end
end

"""
    normalize_kpoint_coordinate(k)

Wrap fractional k-point coordinates into the range [-0.5, 0.5). Accepts a scalar or vector.
"""
function normalize_kpoint_coordinate(x::Real)
    x = x - round(Int, x, RoundNearestTiesUp)
    @assert -0.5 ≤ x < 0.5
    return x
end
normalize_kpoint_coordinate(k::AbstractVector) = normalize_kpoint_coordinate.(k)

function check_kpoint_reduction(
    rotations, time_reversal, is_shift, mesh, grid_address, k_all_reducible, kirreds
)
    all_rotations = time_reversal ? rotations ∪ [SMatrix{3,3}(-I(3))] : rotations
    shift = is_shift .// 2

    for (iks_reducible, k) in zip(k_all_reducible, kirreds), ikred in iks_reducible
        kred = (shift .+ grid_address[ikred]) .// mesh
        found_mapping = any(all_rotations) do W
            # If the difference between kred and W' * k == W^{-1} * k
            # is only integer in fractional reciprocal-space coordinates, then
            # kred and S' * k are equivalent k-points
            all(isinteger, kred - (W * k))
        end
        if !found_mapping
            return false
        end
    end
    return true
end

function get_recip_mesh(atoms::FlexibleSystem, density, mesh)
    if isnothing(mesh)
        real_lattice = ustrip.(stack(cell_vectors(atoms)))
        return get_recip_mesh(real_lattice, density)
    else
        return mesh
    end
end

function get_recip_mesh(real_lattice, density)
    recip_basis_lengths = norm.(eachcol(get_recip_lattice(real_lattice)))
    return ceil.(Int, recip_basis_lengths .* density)
end

"""
    precompute_phases(kpoints, images) -> Matrix{ComplexF64}

Compute Bloch phase factors ``e^{2πi \\mathbf{T} \\cdot \\mathbf{k}}`` for all
image–k-point pairs. Returns a `(nimages, nkpoints)` matrix.
"""
function precompute_phases(ks, ts)
    K = reduce(hcat, ks)
    T = reduce(hcat, ts)

    dots = T' * K

    return cis.(2π .* dots)
end