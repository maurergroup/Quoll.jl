using Configurations
using ..Basis
using ..CoreProjection

@option struct CoreProjectionParams
    core_basis::Vector{BasisMetadata}
    method::CoreProjection.AbstractCoreProjection
end
