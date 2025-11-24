function DeepHOperator(in_operator::BSparseOperator; radii = nothing, hermitian = false, float = Float64)

    in_metadata = get_metadata(in_operator)
    in_atoms = get_atoms(in_operator)
    in_basisset = get_basisset(in_operator)
    in_spins = get_spins(in_operator)
    in_kind = get_kind(in_operator)

    # Convert sparsity
    out_sparsity = convert_sparsity(in_metadata, radii, RealBlockSparsity, hermitian = hermitian)

    # Construct metadata
    out_metadata = DeepHMetadata(in_atoms, out_sparsity, in_basisset, in_spins)
    
    # Convert data
    out_data = convert_operator_data(out_metadata, in_operator, float = float)

    return DeepHOperator(in_kind, out_data, out_metadata)
end

function convert_operator_data(out_metadata::DeepHMetadata, in_operator::BSparseOperator; float = Float64)
    return convert_operator_data(
        get_sparsity(out_metadata),
        get_basisset(out_metadata),
        get_keydata(in_operator),
        get_atoms(in_operator),
        get_sparsity(in_operator),
        DeepHOperator,
        BSparseOperator,
        float,
    )
end

function convert_operator_data(out_sparsity, out_basisset, in_keydata, in_atoms, in_sparsity,
    out_type::Type{DeepHOperator}, in_type::Type{BSparseOperator}, ::Type{T} = Float64) where T

    herm_nonherm_conv = in_sparsity.hermitian && !out_sparsity.hermitian
    species2nbasis = get_species2nbasis(out_basisset)

    shconv = SHConversion(out_type) ∘ inv(SHConversion(in_type))
    orders = precompute_orders(out_basisset, shconv)
    phases = precompute_phases(out_basisset, shconv, DIM = Val(2))

    out_keys = NTuple{5, Int}[]
    out_blocks = Array{T, 2}[]

    for ((iat, jat), out_images_ij) in pairs(convert_to_hermitian(out_sparsity).ij2images)
        zi = out_basisset.atom2species[iat]
        zj = out_basisset.atom2species[jat]

        # Index with chemical species before entering hot loops
        in_keydata_zizj = in_keydata[(iat, jat)]
        orders_zi = orders[zi]
        orders_zj = orders[zj]
        phases_zizj = phases[(zi, zj)]
        
        # Find the mapping between input and output image indices
        in_images_ij = in_sparsity.ij2images[(iat, jat)]
        ilocalout_to_ilocalin_ij = indexin(out_images_ij, in_images_ij)

        for (out_iR, out_R) in enumerate(out_images_ij)
            push!(out_keys, (out_R..., iat, jat))

            block = Array{T, 2}(undef, species2nbasis[zi], species2nbasis[zj])
            in_iR = ilocalout_to_ilocalin_ij[out_iR]
            
            # in_iR can be nothing if in_sparsity is sparser,
            # e.g. obtained from small radii
            if !isnothing(in_iR)
                for jb in axes(block, 2)
                    @inbounds for ib in axes(block, 1)
                        block[ib, jb] = in_keydata_zizj[orders_zi[ib], orders_zj[jb], in_iR] * phases_zizj[ib, jb]
                    end
                end
            end
            push!(out_blocks, block)

        end
    end

    if herm_nonherm_conv
        for ikey in axes(out_keys, 1)
            key = out_keys[ikey]
            out_image, iat, jat = key[1:3], key[4], key[5]

            # On-site RealBlockSparsity is always complete
            iat == jat && continue

            push!(out_keys, (.-out_image..., jat, iat))
            push!(out_blocks, collect(adjoint(out_blocks[ikey])))
        end
    end

    out_data = Dictionary(out_keys, out_blocks)

    return out_data
end
