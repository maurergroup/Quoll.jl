function BSparseOperator(in_operator::FHIaimsCSCOperator; radii = nothing, hermitian = true, float = Float64)

    in_metadata = get_metadata(in_operator)
    in_atoms = get_atoms(in_operator)
    in_basisset = get_basisset(in_operator)
    in_spins = get_spins(in_operator)
    in_kind = get_kind(in_operator)

    # Convert sparsity
    out_sparsity = convert_sparsity(in_metadata, radii, RealBlockSparsity, hermitian = hermitian)

    # Construct metadata
    out_metadata = RealBSparseMetadata(in_atoms, out_sparsity, in_basisset, in_spins)

    # Initialize out_operator with zeros
    out_operator = build_operator(BSparseOperator, in_kind, out_metadata, uninit = false, value = zero(float))

    # Populate out_operator with values from the in_operator
    convert_operator_data!(out_operator, in_operator)

    return out_operator
end

function convert_operator_data!(out_operator::BSparseOperator, in_operator::FHIaimsCSCOperator)

    out_keydata = get_keydata(out_operator)
    out_sparsity = get_sparsity(out_operator)
    out_basisset = get_basisset(out_operator)
    in_data = get_data(in_operator)
    in_sparsity = get_sparsity(in_operator)

    # TODO: We could perform hermitian to hermitian convert_operator_data! and afterwards perform
    # BSparseOperator hermitian -> BSparseOperator non-hermitian conversion.
    # However, one would have to modify BSparseOperator(::FHIaimsCSCOperator)
    # by computing non-hermitian sparsity and metadata after convert_operator_data! and initiating
    # the conversion.
    herm_to_nonherm = in_sparsity.hermitian && !out_sparsity.hermitian
    herm_to_nonherm && throw(error("Hermitian to non-hermitian convert_operator_data! is not implemented"))

    rowval, colcellptr, in_images = in_sparsity.rowval, in_sparsity.colcellptr, in_sparsity.images
    atom2species, out_images = out_basisset.atom2species, out_sparsity.images
    
    atom2offset = get_atom2offset(out_basisset)
    basis2atom = get_basis2atom(out_basisset)
    iexternal2ilocal = get_iexternal2ilocal(in_images, out_sparsity)

    shconv = SHConversion(BSparseOperator) ∘ inv(SHConversion(FHIaimsCSCOperator))

    # Use inv(shconv) for orders and phases because matrix is being converted
    # using `type2` conversion, for more information see tests/unit/shconversion.jl
    orders = precompute_orders(out_basisset, inv(shconv))
    phases = precompute_phases(out_basisset, inv(shconv), DIM = Val(2))

    natoms = length(atom2species)

    # Convert from dictionaries to arrays for faster access
    orders_atomarray = convert_to_atomarray(orders, atom2species)
    phases_atomarray = convert_to_atomarray(phases, atom2species)
    out_keydata_atomarray = convert_to_atomarray(out_keydata, natoms)
    iexternal2ilocal_atomarray = convert_to_atomarray(iexternal2ilocal, natoms)

    for iR in axes(colcellptr, 2)
        R = in_images[iR]
        R ∈ out_images || continue

        for jb in axes(colcellptr, 3)
            jat = basis2atom[jb]
            jb_local = jb - atom2offset[jat]
            orders_jat = orders_atomarray[jat]

            i_index_first = colcellptr[1, iR, jb]
            i_index_last = colcellptr[2, iR, jb]

            @inbounds for i_index in i_index_first:i_index_last
                ib = rowval[i_index]
                iat = basis2atom[ib]
                ib_local = ib - atom2offset[iat]
                orders_iat = orders_atomarray[iat]

                iR_local = iexternal2ilocal_atomarray[iat, jat][iR]
                !isnothing(iR_local) || continue

                out_keydata_atomarray[iat, jat][
                    orders_iat[ib_local],
                    orders_jat[jb_local],
                    iR_local
                ] = in_data[i_index] * phases_atomarray[iat, jat][ib_local, jb_local]
            end
        end
    end

    onsite_ilocal2imlocal = get_onsite_ilocal2imlocal(out_sparsity)
    species2nbasis = get_species2nbasis(out_basisset)

    for iat in 1:natoms
        z = atom2species[iat]
        nbasis = species2nbasis[z]
        out_keydata_iat = out_keydata_atomarray[iat, iat]

        for iR in axes(out_sparsity.ij2images[(iat, iat)], 1)
            imR = onsite_ilocal2imlocal[iat][iR]

            for jb in 1:nbasis
                @inbounds for ib in 1:nbasis
                    jb > ib || continue
                    out_keydata_iat[jb, ib, imR] = out_keydata_iat[ib, jb, iR]
                end
            end

        end
    end

end

function convert_sparsity(in_metadata::FHIaimsCSCMetadata, out_sparsity_type::Type{RealBlockSparsity}; hermitian = true)
    return convert_sparsity(get_sparsity(in_metadata), get_basisset(in_metadata), out_sparsity_type, hermitian = hermitian)
end
