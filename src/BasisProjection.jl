module BasisProjection
export AbstractBasisProjection, IMPLEMENTED_PROJECTIONS, FC99V

abstract type AbstractBasisProjection end

struct FC99V <: AbstractBasisProjection end

const IMPLEMENTED_PROJECTIONS = (FC99V,)
    
end