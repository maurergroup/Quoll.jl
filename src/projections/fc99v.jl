struct FC99V <: AbstractBasisProjection end

basis_projection(::Val{:fc99v}) = FC99V