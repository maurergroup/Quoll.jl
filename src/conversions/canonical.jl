### FOURIER TRANSFORM ###

# There is a method that uses raw data which can be used instead
function fourier_transform_data!(
    ::Type{KDₒᵤₜ}, ::Type{Dₒᵤₜ}, ::Type{KDᵢₙ}, ::Type{Dᵢₙ},
    out_operator::AbstractOperator, in_operator::AbstractOperator,
    phases_k,
) where {
    KDₒᵤₜ<:DataContainer,
    Dₒᵤₜ<:DenseRecipData,
    KDᵢₙ<:CanonicalBlockRealKeyData,
    Dᵢₙ<:CanonicalBlockRealData,
}
    return fourier_transform_data!(Dₒᵤₜ, KDᵢₙ, Dᵢₙ, out_operator, in_operator, phases_k)
end

function fourier_transform_data!(
    ::Type{Dₒᵤₜ}, ::Type{KDᵢₙ}, ::Type{Dᵢₙ},
    out_operator::AbstractOperator, in_operator::AbstractOperator,
    phases_k,
) where {
    Dₒᵤₜ<:DenseRecipData,
    KDᵢₙ<:CanonicalBlockRealKeyData,
    Dᵢₙ<:CanonicalBlockRealData,
}
    out_data = op_data(out_operator)
    out_basisset = op_basisset(out_operator)

    in_keydata = op_keydata(in_operator)
    in_data = op_data(in_operator)
    in_sparsity = op_sparsity(in_operator)

    out_shconv = op_shconv(out_operator)
    in_shconv = op_shconv(in_operator)
    Δshconv = out_shconv ∘ inv(in_shconv)
    shconv_isidentity = Val(isidentity(Δshconv))

    return fourier_transform_data!(
        out_data,
        in_keydata,
        in_data,
        out_basisset,
        in_sparsity,
        Δshconv,
        shconv_isidentity,
        phases_k,
    )
end

function fourier_transform_data!(
    out_data::DenseRecipData{T},
    in_keydata::CanonicalBlockRealKeyData,
    ::CanonicalBlockRealData,
    out_basisset::BasisSetMetadata,
    in_sparsity::BlockRealSparsity,
    Δshconv::SHConvention,
    ::Val{isidentity},
    phases_k,
) where {T,isidentity}
    out_data_body = unwrap_data(out_data)
    in_keydata_body = unwrap_data(in_keydata)

    ilocal2iglobal = get_ilocal2iglobal(in_sparsity)
    atom2offset = get_atom2offset(out_basisset)

    if !isidentity # type2 conversion
        atom2species = out_basisset.atom2species
        shifts = precompute_shifts(out_basisset, inv(Δshconv))
        shphases = precompute_shphases(out_basisset, inv(Δshconv), Val(2))
    end

    for ((iat, jat), in_block) in pairs(in_keydata_body)
        phases_kij = phases_k[ilocal2iglobal[(iat, jat)]]

        ## Method: explicit loops

        if !isidentity
            shifts_zi = shifts[atom2species[iat]]
            shifts_zj = shifts[atom2species[jat]]
            shphases_zizj = shphases[atom2species[iat], atom2species[jat]]
        end

        atom2offset_i = atom2offset[iat]
        atom2offset_j = atom2offset[jat]

        for jb in axes(in_block, 2)
            jb_dense = jb + atom2offset_j
            !isidentity && (jb_dense += shifts_zj[jb])
            for ib in axes(in_block, 1)
                ib_dense = ib + atom2offset_i
                !isidentity && (ib_dense += shifts_zi[ib])

                tmp = zero(T)
                @inbounds for iR in eachindex(phases_kij)
                    tmp = muladd(in_block[ib, jb, iR], phases_kij[iR], tmp)
                end

                out_data_body[ib_dense, jb_dense] = tmp
                !isidentity && (out_data_body[ib_dense, jb_dense] *= shphases_zizj[ib, jb])
            end
        end

        ## Method: matrix-vector product
        # (doesn't include shconv)

        # in_matrix = reshape(in_block, size(in_block, 1) * size(in_block, 2), size(in_block, 3))
        # out_data_body[atom2basis[iat], atom2basis[jat]] = in_matrix * phases_kij

        ## Method: using Tullio
        # (doesn't include shconv)

        # @tullio threads=false grad=false tmp[ib,jb] := begin
        #     @inbounds in_block[ib,jb,R] * phases_kij[R]
        # end
        # out_data_body[atom2basis[iat], atom2basis[jat]] = tmp
    end

    if in_sparsity.hermitian
        for jb_dense in axes(out_data_body, 2), ib_dense in axes(out_data_body, 1)
            if ib_dense < jb_dense
                out_data_body[jb_dense, ib_dense] = conj(out_data_body[ib_dense, jb_dense])
            end
        end
    end
end

### INVERSE FOURIER TRANSFORM ###

function inv_fourier_transform_data!(
    ::Type{KDₒᵤₜ}, ::Type{Dₒᵤₜ}, ::Type{KDᵢₙ}, ::Type{Dᵢₙ},
    out_operator::AbstractOperator, in_operator::AbstractOperator,
    phases_k,
    weight,
) where {
    KDₒᵤₜ<:CanonicalBlockRealKeyData,
    Dₒᵤₜ<:CanonicalBlockRealData,
    KDᵢₙ<:DataContainer,
    Dᵢₙ<:DenseRecipData,
}
    return inv_fourier_transform_data!(
        KDₒᵤₜ, Dₒᵤₜ, Dᵢₙ, out_operator, in_operator, phases_k, weight
    )
end

function inv_fourier_transform_data!(
    ::Type{KDₒᵤₜ}, ::Type{Dₒᵤₜ}, ::Type{Dᵢₙ},
    out_operator::AbstractOperator, in_operator::AbstractOperator,
    phases_k,
    weight,
) where {
    KDₒᵤₜ<:CanonicalBlockRealKeyData,
    Dₒᵤₜ<:CanonicalBlockRealData,
    Dᵢₙ<:DenseRecipData,
}
    out_keydata = op_keydata(out_operator)
    out_data = op_data(out_operator)
    out_sparsity = op_sparsity(out_operator)

    in_basisset = op_basisset(in_operator)
    in_data = op_data(in_operator)

    out_shconv = op_shconv(out_operator)
    in_shconv = op_shconv(in_operator)
    Δshconv = out_shconv ∘ inv(in_shconv)
    shconv_isidentity = Val(isidentity(Δshconv))

    return inv_fourier_transform_data!(
        out_keydata,
        out_data,
        in_data,
        out_sparsity,
        in_basisset,
        Δshconv,
        shconv_isidentity,
        phases_k,
        weight,
    )
end

function inv_fourier_transform_data!(
    out_keydata::CanonicalBlockRealKeyData,
    ::CanonicalBlockRealData,
    in_data::DenseRecipData,
    out_sparsity::BlockRealSparsity,
    in_basisset::BasisSetMetadata,
    Δshconv::SHConvention,
    shconv_isidentity::Val{isidentity},
    phases_k,
    weight,
) where {isidentity}
    out_keydata_body = unwrap_data(out_keydata)
    in_data_body = unwrap_data(in_data)

    ilocal2iglobal = get_ilocal2iglobal(out_sparsity)
    atom2offset = get_atom2offset(in_basisset)

    if !isidentity # type1
        atom2species = in_basisset.atom2species
        shifts = precompute_shifts(in_basisset, Δshconv)
        shphases = precompute_shphases(in_basisset, Δshconv, Val(2))
    end

    for ((iat, jat), out_block) in pairs(out_keydata_body)
        inv_phases_kij = conj(phases_k[ilocal2iglobal[(iat, jat)]])

        ## Method: explicit loops

        atom2offset_i = atom2offset[iat]
        atom2offset_j = atom2offset[jat]

        if !isidentity
            shifts_zi = shifts[atom2species[iat]]
            shifts_zj = shifts[atom2species[jat]]
            shphases_zizj = shphases[atom2species[iat], atom2species[jat]]
        end

        for iR in axes(out_block, 3)
            for jb in axes(out_block, 2)
                jb_dense = jb + atom2offset_j
                !isidentity && (jb_dense += shifts_zj[jb])

                @inbounds for ib in axes(out_block, 1)
                    ib_dense = ib + atom2offset_i
                    !isidentity && (ib_dense += shifts_zi[ib])

                    contribution = (
                        weight *
                        real(in_data_body[ib_dense, jb_dense] * inv_phases_kij[iR])
                    )
                    !isidentity && (contribution *= shphases_zizj[ib, jb])

                    out_block[ib, jb, iR] += contribution
                end
            end
        end

        ## Method: using Tullio

        # in_data_ij = @view(in_data_body[atom2basis[iat], atom2basis[jat]])
        # @tullio threads=false grad=false out_block[ib,jb,iR] += begin
        #     # @inbounds weight * real(in_data_ij[ib,jb] * inv_phases_kij[iR])
        #     weight * real(in_data_ij[ib,jb] * inv_phases_kij[iR])
        # end
    end
end

### CONVERSIONS ###

function convert_operator_data!(
    ::Type{KDₒᵤₜ}, ::Type{Dₒᵤₜ}, ::Type{Dᵢₙ},
    out_operator::AbstractOperator, in_operator::AbstractOperator,
) where {
    KDₒᵤₜ<:CanonicalBlockRealKeyData,
    Dₒᵤₜ<:CanonicalBlockRealData,
    Dᵢₙ<:CSCRealData,
}
    out_hermitian = op_hermicity(out_operator)
    in_hermitian = op_hermicity(in_operator)

    herm_to_herm = in_hermitian && out_hermitian
    if !herm_to_herm
        throw(
            error(
                "only hermitian-to-hermitian conversion is implemented ",
                "for this set of data types",
            ),
        )
    end

    out_keydata = op_keydata(out_operator)
    out_data = op_data(out_operator)
    out_sparsity = op_sparsity(out_operator)
    out_basisset = op_basisset(out_operator)

    in_data = op_data(in_operator)
    in_sparsity = op_sparsity(in_operator)
    in_basisset = op_basisset(in_operator)

    out_shconv = op_shconv(out_operator)
    in_shconv = op_shconv(in_operator)
    Δshconv = out_shconv ∘ inv(in_shconv)
    shconv_isidentity = Val(isidentity(Δshconv))

    return convert_operator_data!(
        out_keydata,
        out_data,
        in_data,
        out_sparsity,
        out_basisset,
        in_sparsity,
        in_basisset,
        Δshconv,
        shconv_isidentity,
    )
end

# TODO:
# - write outer function as for fourier_transform_data
# - make sure this works for differing sparsities
function convert_operator_data!(
    out_keydata::CanonicalBlockRealKeyData,
    ::CanonicalBlockRealData,
    in_data::CSCRealData,
    out_sparsity::BlockRealSparsity,
    out_basisset::BasisSetMetadata,
    in_sparsity::CSCRealSparsity,
    in_basisset::BasisSetMetadata,
    Δshconv::SHConvention,
    ::Val{isidentity},
) where {isidentity}
    out_keydata_body = unwrap_data(out_keydata)
    in_data_body = unwrap_data(in_data)

    rowval = in_sparsity.rowval
    colcellptr = in_sparsity.colcellptr
    in_images = in_sparsity.images
    atom2species = in_basisset.atom2species
    out_images = out_sparsity.images

    atom2offset = get_atom2offset(in_basisset)
    basis2atom = get_basis2atom(in_basisset)
    iexternal2ilocal = get_iexternal2ilocal(in_images, out_sparsity)
    natoms = length(atom2species)

    if !isidentity # type2
        orders = precompute_orders(out_basisset, inv(Δshconv))
        shphases = precompute_shphases(out_basisset, inv(Δshconv), Val(2))
    end

    # Convert from dictionaries to arrays for faster access
    orders_atomarray = convert_to_atomarray(orders, atom2species)
    shphases_atomarray = convert_to_atomarray(shphases, atom2species)
    out_keydata_atomarray = convert_to_atomarray(out_keydata_body, natoms)
    iexternal2ilocal_atomarray = convert_to_atomarray(iexternal2ilocal, natoms)

    for iR in axes(colcellptr, 2)
        R = in_images[iR]

        # out_sparsity could be sparser (e.g. generated with short radii).
        R ∈ out_images || continue

        for jb in axes(colcellptr, 3)
            jat = basis2atom[jb]
            jb_block = jb - atom2offset[jat]
            jb_block2 = isidentity ? jb_block : orders_atomarray[jat][jb_block]

            i_index_first = colcellptr[1, iR, jb]
            i_index_last = colcellptr[2, iR, jb]

            @inbounds for i_index in i_index_first:i_index_last
                ib = rowval[i_index]
                iat = basis2atom[ib]
                ib_block = ib - atom2offset[iat]
                ib_block2 = isidentity ? ib_block : orders_atomarray[iat][ib_block]

                # We had checked whether R is in out_images, whereas here
                # we check whether R is in out_images_ij.
                # This check alone here would be enough, but we can possibly
                # save some additional time if we do the check outside the loop as well
                iR_local = iexternal2ilocal_atomarray[iat, jat][iR]
                !isnothing(iR_local) || continue

                out_keydata_atomarray[iat, jat][ib_block2, jb_block2, iR_local] =
                    in_data_body[i_index]

                if !isidentity
                    out_keydata_atomarray[iat, jat][ib_block2, jb_block2, iR_local] *=
                        shphases_atomarray[iat, jat][ib_block, jb_block]
                end
            end
        end
    end

    onsite_ilocal2imlocal = get_onsite_ilocal2imlocal(out_sparsity)
    species2nbasis = get_species2nbasis(in_basisset)

    # Populating non-hermitian values of on-site blocks. This is because
    # on-site BlockRealSparsity is always non-hermitian (even when hermitian = true).
    # If this was not the case, we would still have to do this at least for
    # R = (0, 0, 0) case, because CSCRealSparsity(hermitian = true) only stores
    # the upper triangle of on-site R = (0, 0, 0) blocks

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
