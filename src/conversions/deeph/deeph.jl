function DeepHOperator(in_operator::RealBSparseOperator; radii = nothing, hermitian = false, T = Float64)

    # Obtain sparsity
    if isnothing(radii)
        if hermitian
            out_sparsity = convert_to_hermitian(in_operator.metadata.sparsity)
        else
            out_sparsity = convert_to_nonhermitian(in_operator.metadata.sparsity)
        end
    else
        @info "Computing sparsity manually from radii"
        out_sparsity = RealBlockSparsity(in_operator.metadata.atoms, radii, hermitian = hermitian)
    end

    # Construct metadata
    out_metadata = DeepHMetadata(
        in_operator.metadata.atoms, out_sparsity, in_operator.metadata.basisset, in_operator.metadata.spins
    )

    # Construct data
    # shconv = SHConversion(DeepHOperator) ∘ inv(SHConversion(RealBSparseOperator))
    # TODO:
    shconv = inv(SHConversion(DeepHOperator) ∘ inv(SHConversion(RealBSparseOperator)))
    shifts, phases = precompute_shiftphases(out_metadata.basisset, shconv)
    
    out_keys = NTuple{5, Int}[]
    out_blocks = Array{T, 2}[]

    Nb_dict = dictionary([
        z => length(basis_species(out_metadata.basisset, z))
        for z in keys(out_metadata.basisset.basis)
    ])

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

    herm_nonherm_conv = in_operator.metadata.sparsity.hermitian && !out_metadata.sparsity.hermitian

    for ((iat, jat), images_ij) in pairs(out_metadata.sparsity.ij2images)
        ij_present = haskey(in_operator.metadata.sparsity.ij2images, (iat, jat))

        zi = out_metadata.basisset.atom2species[iat]
        zj = out_metadata.basisset.atom2species[jat]

        for image in images_ij
            push!(out_keys, (image..., iat, jat))

            block = Array{T, 2}(undef, Nb_dict[zi], Nb_dict[zj])
            if herm_nonherm_conv && (!ij_present || image ∉ in_operator.metadata.sparsity.ij2images[(iat, jat)])
                mt_image = Tuple(-image) # need to check allocations for this, maybe use view
                for jb in axes(block, 2)
                    for ib in axes(block, 1)
                        block[ib + shifts[zi][ib], jb + shifts[zj][jb]] = in_operator[(jat, iat)][
                            jb, ib, Key(mt_image)
                        ] * phases[zi][ib] * phases[zj][jb]
                    end
                end
            else
                t_image = Tuple(image)
                for jb in axes(block, 2)
                    for ib in axes(block, 1)
                        block[ib + shifts[zi][ib], jb + shifts[zj][jb]] = in_operator[(iat, jat)][
                            ib, jb, Key(t_image)
                        ] * phases[zi][ib] * phases[zj][jb]
                    end
                end
            end
            push!(out_blocks, block)

        end
    end
    out_data = Dictionary(out_keys, out_blocks)

    return DeepHOperator(in_operator.kind, out_data, out_metadata)
end