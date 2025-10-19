
# TODO: might need to write a more specialized function which doesn't do conversions for
# each operator separately assuming the metadata is the same
function RealBSparseOperator(in_operator::FHIaimsCSCOperator)
    ## Metadata
    # Convert sparsity
    out_sparsity = RealBlockSparsity(in_operator.metadata.sparsity, in_operator.metadata.basisset)

    # Compute z1z2_ij2interval
    z1z2_ij2interval = compute_z1z2_ij2interval(in_operator.metadata.atoms, out_sparsity)

    out_metadata = RealBSparseMetadata(
        in_operator.metadata.atoms, out_sparsity, in_operator.metadata.basisset, in_operator.metadata.spins, z1z2_ij2interval
    )
    
    ## Data
    # Initialize out_operator with zeros
    out_operator = RealBSparseOperator(in_operator.kind, out_metadata)

    # Populate out_operator with values from the in_operator
    populate!(out_operator, in_operator)

    return out_operator
end

# Probably shouldn't be used directly because this assumes appropriately converted metadata
function populate!(out_operator::RealBSparseOperator, in_operator::FHIaimsCSCOperator)
    return populate!(
        out_operator.keydata,
        out_operator.metadata.sparsity,
        out_operator.metadata.basisset,
        in_operator.data,
        in_operator.metadata.sparsity,
        RealBSparseOperator,
        FHIaimsCSCOperator,
    )
end

# Loop over the CSC sparsity and occupy appropriate values based on block sparsity
function populate!(out_keydata, out_sparsity, out_basisset, in_data, in_sparsity, out_type::Type{RealBSparseOperator}, in_type::Type{FHIaimsCSCOperator})

    atom2basis_offset = get_offsets(out_basisset)
    iglobal2ilocal = get_iglobal2ilocal(out_sparsity)

    # shconv = SHConversion(out_type) ∘ inv(SHConversion(in_type))
    # shifts, phases = precompute_shiftphases(out_basisset, shconv)

    for i_cell in axes(in_sparsity.colcellptr, 2)
        image = in_sparsity.images[i_cell]
        image ∈ out_sparsity.images || continue

        for i_basis_col in axes(in_sparsity.colcellptr, 3)
            i_atom_col = out_basisset.basis2atom[i_basis_col]

            i_index_first = in_sparsity.colcellptr[1, i_cell, i_basis_col]
            i_index_last = in_sparsity.colcellptr[2, i_cell, i_basis_col]

            for i_index in i_index_first:i_index_last
                i_basis_row = in_sparsity.rowval[i_index]
                i_atom_row = out_basisset.basis2atom[i_basis_row]

                i_cell_local = iglobal2ilocal[(i_atom_row, i_atom_col)][i_cell]
                !isnothing(i_cell_local) || continue

                i_basis_row_local = i_basis_row - atom2basis_offset[i_atom_row]
                i_basis_col_local = i_basis_col - atom2basis_offset[i_atom_col]

                out_keydata[(i_atom_row, i_atom_col)][i_basis_row_local, i_basis_col_local, i_cell_local] = in_data[i_index]
            end
        end
    end
end

# struct FHIaimsCSCMetadata{A<:AbstractSystem, E} <: AbstractFHIaimsMetadata
#     atoms::A
#     sparsity::RealCSCSparsity
#     basisset::BasisSetMetadata{E}
#     spinset::Union{SpinsMetadata, Nothing}
#     # TODO: Is Union{SpinsMetadata, Nothing} best approach here?
#     # Making FHIaimsCSCMetadata a parametric type wrt typeof(spins)
#     # or making SpinsMetadata a parametric type might be an overkill
# end

# struct FHIaimsCSCOperator{O<:AbstractOperatorKind, T<:AbstractFloat, A<:AbstractSystem, E} <: AbstractFHIaimsOperator
#     kind::O
#     data::Vector{T}
#     metadata::FHIaimsCSCMetadata{A, E}
# end