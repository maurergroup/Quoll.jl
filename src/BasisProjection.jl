module BasisProjection
export AbstractBasisProjection, PROJ_REGISTRY, FC99V

abstract type AbstractBasisProjection end

struct FC99V <: AbstractBasisProjection end

const PROJ_REGISTRY = (FC99V,)
    
end