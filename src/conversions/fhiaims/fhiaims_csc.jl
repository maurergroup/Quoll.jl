
function RealBSparseOperator(in_operator::FHIaimsCSCOperator; radii = nothing, hermitian = true, float::Type{T} = Float64) where T
# function RealBSparseOperator(in_operator::FHIaimsCSCOperator; radii = nothing, hermitian = true, float = Float64)

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
    # out_operator = RealBSparseOperator(in_kind, out_metadata; float = float)
    out_operator = RealBSparseOperator(in_kind, out_metadata, float)

    # Populate out_operator with values from the in_operator
    populate!(out_operator, in_operator)

    return out_operator
end

# Probably shouldn't be used directly because this assumes appropriately converted metadata
function populate!(out_operator::RealBSparseOperator, in_operator::FHIaimsCSCOperator)
    return populate!(
        get_data(out_operator),
        get_atoms(out_operator),
        get_sparsity(out_operator),
        get_basisset(out_operator),
        get_data(in_operator),
        get_sparsity(in_operator),
        RealBSparseOperator,
        FHIaimsCSCOperator,
    )
end

# Loop over the CSC sparsity and occupy appropriate values based on block sparsity
function populate!(out_data, out_atoms, out_sparsity, out_basisset, in_data, in_sparsity,
    out_type::Type{RealBSparseOperator}, in_type::Type{FHIaimsCSCOperator})
    # TODO: We could perform hermitian to hermitian populate! and afterwards perform
    # RealBSparseOperator hermitian -> RealBSparseOperator non-hermitian conversion.
    # However, one would have to modify RealBSparseOperator(::FHIaimsCSCOperator)
    # by computing non-hermitian sparsity and metadata after populate! and initiating
    # the conversion.
    herm_to_nonherm = in_sparsity.hermitian && !out_sparsity.hermitian
    herm_to_nonherm && throw(error("Hermitian to non-hermitian populate! is not implemented"))
    
    atom2offset = get_atom2offset(out_basisset)
    basis2atom = get_basis2atom(out_basisset)
    iglobal2ilocal = get_iglobal2ilocal(out_sparsity)
    ij2offset = get_ij2offset(out_atoms, out_sparsity)

    shconv = SHConversion(out_type) ∘ inv(SHConversion(in_type))

    # Use inv(shconv) for orders and phases because matrix is being converted
    # using `type2` conversion, for more information see tests/unit/shconversion.jl
    orders = precompute_orders(out_basisset, inv(shconv))
    phases = precompute_phases(out_basisset, inv(shconv), DIM = 2)

    rowval, colcellptr, in_images = in_sparsity.rowval, in_sparsity.colcellptr, in_sparsity.images
    atom2species, out_images = out_basisset.atom2species, out_sparsity.images

    # TODO: The conversion is quite slow mostly because it contains many dict accesses
    # in the loop. I tried to optimise this with more suitable dictionaries below
    # with only limited success. One could instead rewrite those data structures
    # using a direct access table instead (e.g. 118 x 118 chemical elements for z1z2,
    # Nat x Nat for ij) if this or similar conversions become a serious bottleneck.
    # Alternatively they could locally be converted to direct access tables here instead
    # as long as computing them is less computationally expensive than the loop with dicts
    # which is probably always going to be the case.

    # Dicts which will always be small (usually not that many z1z2 pairs)
    out_data_view = Base.ImmutableDict(collect(pairs(out_data))...)
    orders_view = Base.ImmutableDict(collect(pairs(orders))...)
    phases_view = Base.ImmutableDict(collect(pairs(phases))...)

    # Persistent is slower than Dictionaries dict when I'm benchmarking with 2 atoms
    # but in theory it should be faster for larger dicts (more atoms)
    ij2offset_view = Base.PersistentDict(collect(pairs(ij2offset))...)
    iglobal2ilocal_view = Base.PersistentDict(collect(pairs(iglobal2ilocal))...)

    for iR in axes(colcellptr, 2)
        R = in_images[iR]
        R ∈ out_images || continue

        for jb in axes(colcellptr, 3)
            jat = basis2atom[jb]
            zj = atom2species[jat]
            jb_local = jb - atom2offset[jat]

            orders_zj = orders_view[zj]

            i_index_first = colcellptr[1, iR, jb]
            i_index_last = colcellptr[2, iR, jb]

            @inbounds for i_index in i_index_first:i_index_last
                ib = rowval[i_index]

                iat = basis2atom[ib]
                zi = atom2species[iat]
                ib_local = ib - atom2offset[iat]

                orders_zi = orders_view[zi]

                iR_local = iglobal2ilocal_view[(iat, jat)][iR]
                !isnothing(iR_local) || continue
                iR_offset = ij2offset_view[(iat, jat)]

                out_data_view[(zi, zj)][
                    orders_zi[ib_local],
                    orders_zj[jb_local],
                    iR_local + iR_offset
                ] = in_data[i_index] * phases_view[(zi, zj)][ib_local, jb_local]
            end
        end
    end

    n_atoms = length(atom2offset)
    onsite_ilocal2imlocal = get_onsite_ilocal2imlocal(out_sparsity)
    species2nbasis = get_species2nbasis(out_basisset)

    for iat in 1:n_atoms
        z = atom2species[iat]
        nbasis = species2nbasis[z]
        out_data_z = out_data[(z, z)]

        iR_offset = ij2offset_view[(iat, iat)]
        for iR in axes(out_sparsity.ij2images[(iat, iat)], 1)
            imR = onsite_ilocal2imlocal[iat][iR]

            out_iR = iR + iR_offset
            out_imR = imR + iR_offset

            for jb in 1:nbasis
                @inbounds for ib in nbasis
                    jb > ib || continue
                    out_data_z[jb, ib, out_imR] = out_data_z[ib, jb, out_iR]
                end
            end

        end
    end

end
