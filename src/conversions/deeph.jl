function convert_operator_data!(
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

    return convert_operator_data!(
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

function convert_operator_data!(
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

        # herm -> nonherm: (jat, iat) will not be present inside in_sparsity.
        # But generally this is a possibility all the time,
        # e.g. if in_sparsity is generated with very small radii,
        # then there might be some pairs with no images
        (iat, jat) ∉ keys(in_sparsity.ij2images) && continue

        zi = out_basisset.atom2species[iat]
        zj = out_basisset.atom2species[jat]

        # Find the mapping between input and output image indices
        in_images_ij = in_sparsity.ij2images[(iat, jat)]
        ilocalout_to_ilocalin_ij = indexin(out_images_ij, in_images_ij)

        # Index with chemical species before entering hot loops
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
                out_block = out_data_body[key]

                for jb in axes(out_block, 2)
                    jb2 = isidentity ? jb : orders_zj[jb]

                    @inbounds for ib in axes(out_block, 1)
                        ib2 = isidentity ? ib : orders_zi[ib]

                        out_block[ib, jb] = in_keydata_ij[ib2, jb2, in_iR]
                        !isidentity && (out_block[ib, jb] *= shphases_zizj[ib, jb])
                    end
                end
            end
        end
    end

    # In the herm -> nonherm case the lower triangle of out_keydata is not filled
    # because in_sparsity only contains atom pairs belonging to the upper triangle.
    # To alleviate this, the lower triangle of out_keydata is filled here.
    # On-site keys are skipped because on-site blocks in BlockRealSparsity
    # are always non-hermitian (even in the case where hermitian = true)
    herm_nonherm_conv = op_hermicity(in_sparsity) && !op_hermicity(out_sparsity)
    if herm_nonherm_conv
        for key in keys(out_data_body)
            Rx, Ry, Rz, iat, jat = key
            iat == jat && continue

            mkey = (-Rx, -Ry, -Rz, jat, iat)
            out_data_body[mkey] .= adjoint(out_data_body[key])
        end
    end
end
