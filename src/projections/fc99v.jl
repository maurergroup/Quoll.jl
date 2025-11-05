struct FC99V <: AbstractBasisProjection end

get_basis_projection(::Val{:fc99v}) = FC99V