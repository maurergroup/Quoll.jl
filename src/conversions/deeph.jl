function convert_data!(
    ::Type{Dₒᵤₜ}, ::Type{KDᵢₙ}, ::Type{Dᵢₙ},
    out_operator::AbstractOperator, in_operator::AbstractOperator,
) where {
    Dₒᵤₜ<:DeepHBlockRealData,
    KDᵢₙ<:CanonicalBlockRealKeyData,
    Dᵢₙ<:CanonicalBlockRealData,
}
    out_data = op_data(out_operator)
    out_sparsity = op_sparsity(out_operator)
    out_basisset = op_basisset(out_operator)
    out_shconv = op_shconv(out_operator)

    in_keydata = op_keydata(in_operator)
    in_data = op_data(in_operator)
    in_sparsity = op_sparsity(in_operator)
    in_shconv = op_shconv(in_operator)

    Δshconv = out_shconv ∘ inv(in_shconv)
    shconv_isidentity = Val(isidentity(Δshconv))

    return convert_data!(
        out_data,
        in_keydata,
        in_data,
        out_sparsity,
        out_basisset,
        in_sparsity,
        Δshconv,
        shconv_isidentity,
    )
end

function convert_data!(
    out_data::DeepHBlockRealData,
    in_keydata::CanonicalBlockRealKeyData,
    ::CanonicalBlockRealData,
    out_sparsity::BlockRealSparsity,
    out_basisset::BasisSetMetadata,
    in_sparsity::BlockRealSparsity,
    Δshconv::SHConvention,
    ::Val{isidentity},
) where {isidentity}
    out_data_body = unwrap_data(out_data)
    in_keydata_body = unwrap_data(in_keydata)

    if !isidentity # type1
        orders = precompute_orders(out_basisset, Δshconv)
        shphases = precompute_shphases(out_basisset, Δshconv, Val(2))
    end

    for ((iat, jat), out_images_ij) in pairs(out_sparsity.ij2images)

        # herm -> nonherm: (jat, iat) might not be present inside in_sparsity.
        # But generally this is a possibility all the time,
        # e.g. if in_sparsity is generated with very small radii,
        # then there might be some pairs with no images
        !haskey(in_sparsity.ij2images, (iat, jat)) && continue

        zi = out_basisset.atom2species[iat]
        zj = out_basisset.atom2species[jat]

        # Find the mapping between input and output image indices
        in_images_ij = in_sparsity.ij2images[(iat, jat)]
        ilocalout_to_ilocalin_ij = indexin(out_images_ij, in_images_ij)

        # Index before entering hot loops
        in_keydata_ij = in_keydata_body[(iat, jat)]

        if !isidentity
            orders_zi = orders[zi]
            orders_zj = orders[zj]
            shphases_zizj = shphases[(zi, zj)]
        end

        for (out_iR, out_R) in enumerate(out_images_ij)
            in_iR = ilocalout_to_ilocalin_ij[out_iR]

            # in_iR can be nothing if in_sparsity is sparser,
            # e.g. obtained from small radii
            if !isnothing(in_iR)
                key = (out_R..., iat, jat)
                out_data_ij = out_data_body[key]

                for jb in axes(out_data_ij, 2)
                    jb2 = isidentity ? jb : orders_zj[jb]

                    @inbounds for ib in axes(out_data_ij, 1)
                        ib2 = isidentity ? ib : orders_zi[ib]

                        out_data_ij[ib, jb] = in_keydata_ij[ib2, jb2, in_iR]
                        !isidentity && (out_data_ij[ib, jb] *= shphases_zizj[ib, jb])
                    end
                end
            end
        end
    end

    # In the herm -> nonherm case the lower triangle ij pairs of out_data
    # are not filled because in_sparsity only contains atom pairs belonging
    # to the upper triangle ij pairs.
    # To alleviate this, the lower triangle ij pairs of out_data are filled here.
    # On-site keys are skipped because on-site blocks in BlockRealSparsity
    # are always non-hermitian (even in the case where hermitian = true)
    herm_nonherm_conv = op_hermicity(in_sparsity) && !op_hermicity(out_sparsity)
    if herm_nonherm_conv
        # TODO: write as an external method
        for key in keys(out_data_body)
            Rx, Ry, Rz, iat, jat = key
            iat == jat && continue

            mkey = (-Rx, -Ry, -Rz, jat, iat)
            out_data_body[mkey] .= adjoint(out_data_body[key])
        end
    end
end

function convert_data!(
    ::Type{KDₒᵤₜ}, ::Type{Dₒᵤₜ}, ::Type{Dᵢₙ},
    out_operator::AbstractOperator, in_operator::AbstractOperator,
) where {
    KDₒᵤₜ<:CanonicalBlockRealKeyData,
    Dₒᵤₜ<:CanonicalBlockRealData,
    Dᵢₙ<:DeepHBlockRealData,
}
    out_keydata = op_keydata(out_operator)
    out_data = op_data(out_operator)
    out_sparsity = op_sparsity(out_operator)
    out_basisset = op_basisset(out_operator)
    out_shconv = op_shconv(out_operator)

    in_data = op_data(in_operator)
    in_sparsity = op_sparsity(in_operator)
    in_shconv = op_shconv(in_operator)

    Δshconv = out_shconv ∘ inv(in_shconv)
    shconv_isidentity = Val(isidentity(Δshconv))

    return convert_data!(
        out_keydata,
        out_data,
        in_data,
        out_sparsity,
        out_basisset,
        in_sparsity,
        Δshconv,
        shconv_isidentity,
    )
end

function convert_data!(
    out_keydata::CanonicalBlockRealKeyData,
    ::CanonicalBlockRealData,
    in_data::DeepHBlockRealData,
    out_sparsity::BlockRealSparsity,
    out_basisset::BasisSetMetadata,
    in_sparsity::BlockRealSparsity,
    Δshconv::SHConvention,
    ::Val{isidentity},
) where {isidentity}
    out_keydata_body = unwrap_data(out_keydata)
    in_data_body = unwrap_data(in_data)

    if !isidentity # type1
        orders = precompute_orders(out_basisset, Δshconv)
        shphases = precompute_shphases(out_basisset, Δshconv, Val(2))
    end

    for ((iat, jat), out_images_ij) in pairs(out_sparsity.ij2images)

        # herm -> nonherm: (jat, iat) might not be present inside in_sparsity.
        # But generally this is a possibility all the time,
        # e.g. if in_sparsity is generated with very small radii,
        # then there might be some pairs with no images
        !haskey(in_sparsity.ij2images, (iat, jat)) && continue

        zi = out_basisset.atom2species[iat]
        zj = out_basisset.atom2species[jat]

        # Index before entering hot loops
        out_keydata_ij = out_keydata_body[(iat, jat)]

        if !isidentity
            orders_zi = orders[zi]
            orders_zj = orders[zj]
            shphases_zizj = shphases[(zi, zj)]
        end

        for (out_iR, out_R) in enumerate(out_images_ij)
            key = (out_R..., iat, jat)

            # key might not be present if in_sparsity is sparser,
            # e.g. obtained from small radii
            if haskey(in_data_body, key)
                in_data_ij = in_data_body[key]

                for jb in axes(out_keydata_ij, 2)
                    jb2 = isidentity ? jb : orders_zj[jb]

                    @inbounds for ib in axes(out_keydata_ij, 1)
                        ib2 = isidentity ? ib : orders_zi[ib]

                        out_keydata_ij[ib, jb, out_iR] = in_data_ij[ib2, jb2]
                        if !isidentity
                            out_keydata_ij[ib, jb, out_iR] *= shphases_zizj[ib, jb]
                        end
                    end
                end
            end
        end
    end

    # In the herm -> nonherm case the lower triangle ij pairs of out_keydata
    # are not filled because in_sparsity only contains atom pairs belonging
    # to the upper triangle ij pairs.
    # To alleviate this, the lower triangle ij pairs of out_keydata are filled here.
    # On-site keys are skipped because on-site blocks in BlockRealSparsity
    # are always non-hermitian (even in the case where hermitian = true)
    herm_nonherm_conv = op_hermicity(in_sparsity) && !op_hermicity(out_sparsity)
    if herm_nonherm_conv
        ilocal2imlocal = get_ilocal2imlocal(out_sparsity)
        for ((iat, jat), out_images_ij) in pairs(out_sparsity.ij2images)
            # Loop only over 'U' pairs
            iat < jat || continue

            ilocal2imlocal_ij = ilocal2imlocal[(iat, jat)]
            out_keydata_ij = out_keydata_body[(iat, jat)]
            out_keydata_ji = out_keydata_body[(jat, iat)]

            for iR in axes(out_images_ij, 1)
                imR = ilocal2imlocal_ij[iR]
                for jb in axes(out_keydata_ij)
                    @inbounds for ib in axes(out_keydata_ij)
                        out_keydata_ji[jb, ib, imR] = conj(out_keydata_ij[ib, jb, iR])
                    end
                end
            end
        end
    end
end

function convert_spins_source(
    in_spins::SpinsMetadata, ::BasisSetMetadata, ::DeepHSource, ::CanonicalSource
)
    return in_spins
end