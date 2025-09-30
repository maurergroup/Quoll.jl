module BasisProjection
export AbstractBasisProjection, FC99V, implemented_basis_projection

abstract type AbstractBasisProjection end

struct FC99V <: AbstractBasisProjection end

function implemented_basis_projection end

implemented_basis_projection(::Val{:fc99v}) = FC99V
    
end