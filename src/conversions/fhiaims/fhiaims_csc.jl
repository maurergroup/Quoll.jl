
function RealBSparseOperator(in_operator::FHIaimsCSCOperator)
    ## Metadata
    # Convert sparsity
    out_sparsity = RealBlockSparsity(in_operator.sparsity, in_operator.basis)

    # Convert basis set IF we decide to store
    # spherical harmonics convention inside it

    # Compute z1z2_ij2interval
    z1z2_ij2interval = compute_z1z2_ij2interval(in_operator.atoms, out_sparsity)

    metadata = RealBSparseMetadata(
        in_operator.atoms, out_sparsity, in_operator.basisset, in_operator.spinset, z1z2_ij2interval
    )
    
    ## Data
    # Initialize bsparse_operator with zeros

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