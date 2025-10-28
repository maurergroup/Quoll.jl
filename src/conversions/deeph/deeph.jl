function DeepHOperator(in_operator::RealBSparseOperator; radii = nothing, hermitian = false, T = Float64)

    # in_sparsity = get_sparsity(in_operator)
    in_metadata = get_metadata(in_operator)
    in_atoms = get_atoms(in_operator)
    in_basisset = get_basisset(in_operator)
    in_spins = get_spins(in_operator)
    in_kind = get_kind(in_operator)

    # Convert sparsity
    # out_sparsity = convert_sparsity(in_sparsity, DeepHOperator, in_atoms = in_atoms, radii = radii, hermitian = hermitian)
    out_sparsity = convert_sparsity(in_metadata, radii, RealBlockSparsity; hermitian = hermitian)

    # Construct metadata
    out_metadata = DeepHMetadata(
        in_atoms, out_sparsity, in_basisset, in_spins
    )
    
    # Convert data
    out_data = convert_operator_data(out_metadata, in_operator, T = T)

    return DeepHOperator(in_kind, out_data, out_metadata)
end

function convert_operator_data(out_metadata::DeepHMetadata, in_operator::RealBSparseOperator; T = Float64)
    return convert_operator_data(
        get_sparsity(out_metadata),
        get_basisset(out_metadata),
        get_keydata(in_operator),
        get_sparsity(in_operator),
        DeepHOperator,
        RealBSparseOperator
    )
end

function convert_operator_data(out_sparsity, out_basisset, in_keydata, in_sparsity, out_type::Type{DeepHOperator}, in_type::Type{RealBSparseOperator}; T = Float64)

    # hermitian -> non-hermitian:
    # - Can be checked if that's the case via hermitian fields in sparsities
    # - Type 1 hermitian sparsity
    #   - ij: [R1, R2]
    #    <ji: ...>: not present
    #   - Would only require to check whether ij is present, if not then (ji, -R)
    # - Type 2 hermitian sparsity
    #   - ij: [R1]
    #     ji: [-R2]
    #   - Even if ij is present, any pair of atoms might not contain all the images e.g. (ij, R2),
    #     in which case it would be required to obtain the block by taking the transpose of the image block
    #     e.g. via (ji, -R2)

    herm_nonherm_conv = in_sparsity.hermitian && !out_sparsity.hermitian

    Nb_dict = dictionary([
        z => length(basis_species(out_basisset, z))
        for z in keys(out_basisset.basis)
    ])

    shconv = SHConversion(out_type) ∘ inv(SHConversion(in_type))
    shifts, phases = precompute_shiftphases(out_basisset, shconv)

    out_keys = NTuple{5, Int}[]
    out_blocks = Array{T, 2}[]

    for ((iat, jat), images_ij) in pairs(out_sparsity.ij2images)
        ij_present = haskey(in_sparsity.ij2images, (iat, jat))

        zi = out_basisset.atom2species[iat]
        zj = out_basisset.atom2species[jat]

        for image in images_ij
            push!(out_keys, (image..., iat, jat))

            block = Array{T, 2}(undef, Nb_dict[zi], Nb_dict[zj])
            if herm_nonherm_conv && (!ij_present || image ∉ in_sparsity.ij2images[(iat, jat)])
                mt_image = Tuple(-image) # need to check allocations for this, maybe use view
                for jb in axes(block, 2)
                    @inbounds for ib in axes(block, 1)
                        block[ib, jb] = in_keydata[(jat, iat)][
                            jb + shifts[zj][jb],
                            ib + shifts[zi][ib],
                            Key(mt_image)
                        ] * phases[zi][ib] * phases[zj][jb]
                    end
                end
            else
                t_image = Tuple(image)
                for jb in axes(block, 2)
                    @inbounds for ib in axes(block, 1)
                        block[ib, jb] = in_keydata[(iat, jat)][
                            ib + shifts[zi][ib],
                            jb + shifts[zj][jb],
                            Key(t_image)
                        ] * phases[zi][ib] * phases[zj][jb]
                    end
                end
            end
            push!(out_blocks, block)

        end
    end
    out_data = Dictionary(out_keys, out_blocks)

    return out_data
end
