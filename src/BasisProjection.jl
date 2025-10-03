module BasisProjection
export AbstractBasisProjection, FC99V, basis_projection

abstract type AbstractBasisProjection end

struct FC99V <: AbstractBasisProjection end

function basis_projection end

basis_projection(::Val{:fc99v}) = FC99V
    
end