
# function DeepHOperator(in_operator::RealBSparseOperator; radii = nothing, hermitian = false)

#     # Obtain sparsity
#     if isnothing(radii)
#         if hermitian
#             out_sparsity = convert_to_hermitian(in_operator.sparsity)
#         else
#             out_sparsity = convert_to_nonhermitian(in_operator.sparsity)
#         end
#     else
#         @info "Computing sparsity manually from radii"
#         out_sparsity = RealBlockSparsity(in_operator.metadata.atoms, radii, hermitian = hermitian)
#     end

#     # Construct metadata
#     out_metadata = DeepHMetadata(
#         in_operator.metadata.atoms, out_sparsity, in_operator.metadata.basisset, in_operator.metadata.spins
#     )

#     # Construct data
    

# end