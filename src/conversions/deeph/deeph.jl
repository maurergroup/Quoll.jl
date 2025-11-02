function DeepHOperator(in_operator::RealBSparseOperator; radii = nothing, hermitian = false, float = Float64)

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

# TODO: This is a bit redundant now
function convert_operator_data(out_metadata::DeepHMetadata, in_operator::RealBSparseOperator; float = Float64)
    return convert_operator_data(
        get_sparsity(out_metadata),
        get_basisset(out_metadata),
        get_data(in_operator),
        get_atoms(in_operator),
        get_sparsity(in_operator),
        DeepHOperator,
        RealBSparseOperator,
        float = float,
    )
end

# Might be faster to make it hermitian first in the case of herm_nonherm_conv
# and then simply add additional blocks
function convert_operator_data(out_sparsity, out_basisset, in_data, in_atoms, in_sparsity,
    out_type::Type{DeepHOperator}, in_type::Type{RealBSparseOperator}; float::Type{T} = Float64) where T

    herm_nonherm_conv = in_sparsity.hermitian && !out_sparsity.hermitian

    species2nbasis = get_species2nbasis(out_basisset)
    z1z2_ij2offset = get_z1z2_ij2offset(in_atoms, in_sparsity) # This is zeros for SiC, I think I need a harder example for regression tests

    shconv = SHConversion(out_type) ∘ inv(SHConversion(in_type))
    orders = precompute_orders(out_basisset, shconv)
    phases = precompute_phases(out_basisset, shconv, DIM = 2)

    out_keys = NTuple{5, Int}[]
    out_blocks = Array{T, 2}[]

    # Converting to hermitian locally because then it works for all cases
    # - hermitian -> hermitian: the types of hermicities are different and then it's gg
    # This would be alleviated if we normalize hermicity
    # But if we normalize hermicity that means we need to change the operator itself as well
    # But that might work ok for RealBSparseOperator
    # "@assert convert_to_hermitian(out_sparsity) == in_sparsity"
    # TODO:
    # Case: we assume i ≤ j in hermicity which works out of the box if we convert from FHIaimsCSCOperator
    #       but there might be conversions where the RealBSparsity is different. What then?
    #       We can either:
    #       - Map values to ji-R during conversion if i ≥ j (more efficient)
    #       - Provide a method which would convert the operator into correct sparsity type
    # Also I did not consider the case when cutoff radii don't match
    # That would be easy to alleviate because we would simply fill the block with zeros
    # if that ijR is not present in the in_sparsity

    for ((iat, jat), out_images_ij) in pairs(convert_to_hermitian(out_sparsity).ij2images)
        zi = out_basisset.atom2species[iat]
        zj = out_basisset.atom2species[jat]

        # Index with chemical species before entering hot loops
        in_data_zizj = in_data[(zi, zj)]
        orders_zi = orders[zi]
        orders_zj = orders[zj]
        phases_zizj = phases[(zi, zj)]
        
        # Find the mapping between input and output image indices
        in_images_ij = in_sparsity.ij2images[(iat, jat)]
        ilocalout_to_ilocalin_ij::Vector{Int} = indexin(out_images_ij, in_images_ij)

        # Find the offset for i_image caused by other ij pairs in the same array
        R_offset = z1z2_ij2offset[(zi, zj)][(iat, jat)]

        for (out_i_image, out_image) in enumerate(out_images_ij)
            push!(out_keys, (out_image..., iat, jat))

            block = Array{T, 2}(undef, species2nbasis[zi], species2nbasis[zj])
            i_image = ilocalout_to_ilocalin_ij[out_i_image] + R_offset
            
            for jb in axes(block, 2)
                @inbounds for ib in axes(block, 1)
                    block[ib, jb] = in_data_zizj[orders_zi[ib], orders_zj[jb], i_image] * phases_zizj[ib, jb]
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
