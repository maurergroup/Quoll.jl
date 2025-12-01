### FOURIER TRANSFORMS ###

# TODO: the same function could be used for both DenseOperator and FHIaimsDenseOperator
# (except for the data I would need to use SH conversion)
# TODO: modify to enable float conversions
function fourier_transform(in_operator::BSparseOperator, kpoints, phases, ::Type{DenseOperator}; float = Float64)
    in_metadata = get_metadata(in_operator)
    in_sparsity = get_sparsity(in_operator)
    in_atoms = get_atoms(in_operator)
    in_basisset = get_basisset(in_operator)
    in_spins = get_spins(in_operator)
    in_kind = get_kind(in_operator)

    # Convert metadata
    out_sparsity = convert_sparsity(in_metadata, RecipDenseSparsity, hermitian = in_sparsity.hermitian)
    out_metadata_list = [
        RecipDenseMetadata(in_atoms, out_sparsity, in_basisset, in_spins, SVector{3}(kpoint))
        for kpoint in kpoints
    ]
    
    # Initialize out_operator with zeros
    out_operator_list = [
        build_operator(DenseOperator, in_kind, out_metadata; uninit = false, value = zero(Complex{float}))
        for out_metadata in out_metadata_list
    ]

    # Perform Fourier transform
    for (out_operator, phases_k) in zip(out_operator_list, phases)
        fourier_transform_data!(out_operator, in_operator, phases_k)
    end
    
    return length(out_operator_list) == 1 ? first(out_operator_list) : out_operator_list
end

# TODO: write another method but with shifts (or alternatively write the method with shifts
# but force the compiler to remove shifts for matching SH convention)
function fourier_transform_data!(out_operator::DenseOperator, in_operator::BSparseOperator, phases_k)
    in_keydata = get_keydata(in_operator)
    in_sparsity = get_sparsity(in_operator)
    out_data = get_data(out_operator)
    out_basisset = get_basisset(out_operator)
    out_float = get_float(out_operator)

    ilocal2iglobal = get_ilocal2iglobal(in_sparsity)
    atom2offset = get_atom2offset(out_basisset)

    for ((iat, jat), in_block) in pairs(in_keydata)

        phases_kij = phases_k[ilocal2iglobal[(iat, jat)]]
        atom2offset_i = atom2offset[iat]
        atom2offset_j = atom2offset[jat]

        for jb in axes(in_block, 2)
            jb_dense = jb + atom2offset_j
            for ib in axes(in_block, 1)
                ib_dense = ib + atom2offset_i

                tmp = zero(out_float)
                @inbounds for iR in axes(in_block, 3)
                    tmp += in_block[ib, jb, iR] * phases_kij[iR]
                end

                out_data[ib_dense, jb_dense] = tmp
            end
        end

        # @tullio tmp[ib,jb] := in_block[ib,jb,R] * phases_kij[R]
        # out_data[atom2basis[iat], atom2basis[jat]] = tmp
    end

    # TODO: I could use Hermitian array instead, probably directly in build stage
    if in_sparsity.hermitian
        for jb_dense in 1:size(out_data, 1)
            for ib_dense in 1:jb_dense
                out_data[jb_dense, ib_dense] = conj(out_data[ib_dense, jb_dense])
            end
        end
    end
end

# Only a contribution to out_operator because in_operator contains only a single k-point
function inv_fourier_transform_data!(out_operator::BSparseOperator, in_operator::DenseOperator, phases_k, weight)
    out_basisset = get_basisset(out_operator)
    out_keydata = get_keydata(out_operator)
    out_sparsity = get_sparsity(out_operator)
    in_data = get_data(in_operator)

    ilocal2iglobal = get_ilocal2iglobal(out_sparsity)
    atom2offset = get_atom2offset(out_basisset)

    for ((iat, jat), out_block) in pairs(out_keydata)

        # Assuming phases is a vector here (i.e. for a single k-point)
        inv_phases_kij = conj(phases_k[ilocal2iglobal[(iat, jat)]])
        atom2offset_i = atom2offset[iat]
        atom2offset_j = atom2offset[jat]
        
        for iR in axes(out_block, 3)
            for jb in axes(out_block, 2)
                jb_dense = jb + atom2offset_j
                for ib in axes(out_block, 1)
                    ib_dense = ib + atom2offset_i
                    out_block[ib, jb, iR] += weight * real(in_data[ib_dense, jb_dense] * inv_phases_kij[iR])
                end
            end
        end
    end

end