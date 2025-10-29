
function RealBSparseOperator(in_operator::FHIaimsCSCOperator; radii = nothing, hermitian = true, float = Float64)

    in_metadata = get_metadata(in_operator)
    in_atoms = get_atoms(in_operator)
    in_basisset = get_basisset(in_operator)
    in_spins = get_spins(in_operator)
    in_kind = get_kind(in_operator)

    # Convert sparsity
    out_sparsity = convert_sparsity(in_metadata, radii, RealBlockSparsity, hermitian = hermitian)

    # Compute z1z2_ij2interval
    z1z2_ij2interval = get_z1z2_ij2interval(in_atoms, out_sparsity)

    # Construct metadata
    out_metadata = RealBSparseMetadata(in_atoms, out_sparsity, in_basisset, in_spins, z1z2_ij2interval)

    # Initialize out_operator with zeros
    out_operator = RealBSparseOperator(in_kind, out_metadata; float = float)

    # Populate out_operator with values from the in_operator
    populate!(out_operator, in_operator)

    return out_operator
end

# Probably shouldn't be used directly because this assumes appropriately converted metadata
function populate!(out_operator::RealBSparseOperator, in_operator::FHIaimsCSCOperator)
    return populate!(
        get_keydata(out_operator),
        get_sparsity(out_operator),
        get_basisset(out_operator),
        get_data(in_operator),
        get_sparsity(in_operator),
        RealBSparseOperator,
        FHIaimsCSCOperator,
    )
end

# Loop over the CSC sparsity and occupy appropriate values based on block sparsity
function populate!(out_keydata, out_sparsity, out_basisset, in_data, in_sparsity,
    out_type::Type{RealBSparseOperator}, in_type::Type{FHIaimsCSCOperator})
    # TODO: We could perform hermitian to hermitian populate! and afterwards perform
    # RealBSparseOperator hermitian -> RealBSparseOperator non-hermitian conversion.
    # However, one would have to modify RealBSparseOperator(::FHIaimsCSCOperator)
    # by computing non-hermitian sparsity and metadata after populate! and initiating
    # the conversion.
    herm_to_nonherm = in_sparsity.hermitian && !out_sparsity.hermitian
    herm_to_nonherm && throw(error("Hermitian to non-hermitian populate! is not implemented"))
    
    atom2basis_offset = get_offsets(out_basisset)
    basis2atom = get_basis2atom(out_basisset)
    iglobal2ilocal = get_iglobal2ilocal(out_sparsity)
    shconv = SHConversion(out_type) ∘ inv(SHConversion(in_type))

    # Use inv(shconv) for shifts and phases because matrix is being reordered
    # using `type2` reordering, for more information see tests/unit/shconversion.jl
    shifts, phases = precompute_shiftphases(out_basisset, inv(shconv))

    rowval, colcellptr, in_images = in_sparsity.rowval, in_sparsity.colcellptr, in_sparsity.images
    atom2species, out_images = out_basisset.atom2species, out_sparsity.images

    for i_cell in axes(colcellptr, 2)
        image = in_images[i_cell]
        image ∈ out_images || continue

        for i_basis_col in axes(colcellptr, 3)
            i_atom_col = basis2atom[i_basis_col]

            i_index_first = colcellptr[1, i_cell, i_basis_col]
            i_index_last = colcellptr[2, i_cell, i_basis_col]

            @inbounds for i_index in i_index_first:i_index_last
                i_basis_row = rowval[i_index]
                i_atom_row = basis2atom[i_basis_row]

                i_cell_local = iglobal2ilocal[(i_atom_row, i_atom_col)][i_cell]
                !isnothing(i_cell_local) || continue

                i_basis_row_local = i_basis_row - atom2basis_offset[i_atom_row]
                i_basis_col_local = i_basis_col - atom2basis_offset[i_atom_col]

                out_keydata[(i_atom_row, i_atom_col)][
                    i_basis_row_local + shifts[atom2species[i_atom_row]][i_basis_row_local],
                    i_basis_col_local + shifts[atom2species[i_atom_col]][i_basis_col_local],
                    i_cell_local
                ] = (
                    in_data[i_index]
                    * phases[atom2species[i_atom_row]][i_basis_row_local]
                    * phases[atom2species[i_atom_col]][i_basis_col_local]
                )
            end
        end
    end

    # Fill in the 'L' part of on-site blocks which, due of Hermitian CSC sparsity,
    # were not filled in the loop above
    n_atoms = length(atom2basis_offset)
    onsite_ilocal2milocal = get_onsite_ilocal2milocal(out_sparsity)
    for iat in 1:n_atoms
        for iR in axes(out_keydata[(iat, iat)], 3)
            miR = onsite_ilocal2milocal[iat][iR]

            for iν in axes(out_keydata[(iat, iat)], 2)
                @inbounds for iμ in axes(out_keydata[(iat, iat)], 1)
                    iν > iμ || continue
                    out_keydata[(iat, iat)][iν, iμ, miR] = out_keydata[(iat, iat)][iμ, iν, iR]
                end
            end

        end
    end

end
