
function RealBSparseOperator(in_operator::FHIaimsCSCOperator; radii = nothing, hermitian = true, float = Float64)

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

    rowval, colcellptr, in_images = in_sparsity.rowval, in_sparsity.colcellptr, in_sparsity.images
    atom2species, out_images = out_basisset.atom2species, out_sparsity.images
    
    atom2offset = get_atom2offset(out_basisset)
    basis2atom = get_basis2atom(out_basisset)
    iglobal2ilocal = get_iglobal2ilocal(out_sparsity)

    shconv = SHConversion(out_type) ∘ inv(SHConversion(in_type))

    # Use inv(shconv) for orders and phases because matrix is being converted
    # using `type2` conversion, for more information see tests/unit/shconversion.jl
    orders = precompute_orders(out_basisset, inv(shconv))
    phases = precompute_phases(out_basisset, inv(shconv), DIM = Val(2))

    natoms = length(atom2species)

    # Convert from dictionaries to arrays for faster access
    orders_atomarray = convert_to_atomarray(orders, atom2species)
    phases_atomarray = convert_to_atomarray(phases, atom2species)
    out_keydata_atomarray = convert_to_atomarray(out_keydata, natoms)
    iglobal2ilocal_atomarray = convert_to_atomarray(iglobal2ilocal, natoms)

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

                iR_local = iglobal2ilocal_atomarray[iat, jat][iR]
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
                @inbounds for ib in nbasis
                    jb > ib || continue
                    out_keydata_iat[jb, ib, imR] = out_data_iat[ib, jb, iR]
                end
            end

        end
    end

end

function convert_sparsity(in_metadata::FHIaimsCSCMetadata, out_sparsity_type::Type{RealBlockSparsity}; hermitian = true)
    return convert_sparsity(get_sparsity(in_metadata), get_basisset(in_metadata), out_sparsity_type, hermitian = hermitian)
end
