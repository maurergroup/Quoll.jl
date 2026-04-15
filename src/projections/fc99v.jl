"""
    FC99V <: AbstractBasisProjection

Frozen-core method with a valence correction [1]. Similar to [`LaikovCore`](@ref) but uses
the raw core-valence overlap `S₁₂` directly instead of the orthogonalised `Ŝ₁₂ = S₁₁⁻¹ S₁₂`:

    Ō₂₂ = O₂₂ + A₂₂ + A₂₂'

where `A₂₂ = S₁₂' * (½ O₁₁ S₁₂ - O₁₂)`.

Requires exactly one `Overlap` operator among the input operators.

# References
[1] doi.org/10.1063/5.0050296
"""
struct FC99V <: AbstractBasisProjection end

function compute_valence_data(
    operators, core_mask_list, valence_mask_list, method::FC99V
)
    # Check that all operators are hermitian
    @argcheck all(operator -> op_hermicity(operator), operators)

    # Requires exactly one overlap matrix for computing the projection
    # (if more than two than it's ambiguous which one to use)
    i_overlap = findall(operator -> op_kind(operator) isa Overlap, operators)
    @argcheck length(i_overlap) == 1

    overlap = operators[first(i_overlap)]
    core_mask_overlap = core_mask_list[first(i_overlap)]
    valence_mask_overlap = valence_mask_list[first(i_overlap)]

    M_overlap = typeof(op_metadata(overlap))
    S₁₂ = unwrap_data(op_data(overlap))[core_mask_overlap, valence_mask_overlap]
    overlap_data₁₂ = wrap_data(M_overlap, S₁₂)

    return [
        compute_valence_data(
            typeof(op_metadata(operator)),
            op_data(operator),
            overlap_data₁₂,
            core_mask,
            valence_mask,
            method,
        )
        for (operator, core_mask, valence_mask) in
        zip(operators, core_mask_list, valence_mask_list)
    ]
end

function compute_valence_data(
    ::Type{M},
    operator_data::DenseRecipData{T},
    overlap_data₁₂::DenseRecipData{T},
    core_mask::BitVector,
    valence_mask::BitVector,
    ::FC99V,
) where {M<:AbstractMetadata,T}
    S₁₂ = unwrap_data(overlap_data₁₂)

    O₁₁ = unwrap_data(operator_data)[core_mask, core_mask]
    O₁₂ = unwrap_data(operator_data)[core_mask, valence_mask]
    O₂₂ = unwrap_data(operator_data)[valence_mask, valence_mask]

    α = T(0.5)
    β = T(-1.0)
    LinearAlgebra.BLAS.hemm!('L', 'U', α, O₁₁, S₁₂, β, O₁₂)

    # A₂₂ = S₁₂' * (0.5 * O₁₁ * S₁₂ - O₁₂)
    A₂₂ = LinearAlgebra.BLAS.gemm('C', 'N', S₁₂, O₁₂)

    Ō₂₂ = O₂₂ + A₂₂ + A₂₂'

    return wrap_data(M, Ō₂₂)
end

### PARSING ###

get_basis_projection(::Val{:fc99v}) = FC99V