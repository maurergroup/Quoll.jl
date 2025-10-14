
function RealBSparseOperator(csc_operator::FHIaimsCSCOperator)
    ## Metadata
    # Convert sparsity
    block_sparsity = RealBlockSparsity(csc_operator.sparsity, csc_operator.basis)

    # Convert basis set IF we decide to store
    # spherical harmonics convention inside it
    
    ## Data
    # Initialize bsparse_operator with zeros
    
    # 

end