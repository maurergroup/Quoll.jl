"""
    LaikovCore <: AbstractBasisProjection

Core-projection method from Ref [1]. An application of this method can be found in Ref [2].
Computes the valence-only operator block by subtracting core contributions using the inverse
of the core-core overlap `S₁₁⁻¹`:

    Ō₂₂ = O₂₂ + A₂₂ + A₂₂'

where `A₂₂ = Ŝ₁₂' * (½ O₁₁ Ŝ₁₂ - O₁₂)` and `Ŝ₁₂ = S₁₁⁻¹ S₁₂`.

Requires exactly one `Overlap` operator among the input operators.

# References
[1] doi.org/10.1002/qua.22767
[2] doi.org/10.1038/s41524-026-02020-1
"""
struct LaikovCore <: AbstractBasisProjection end

"""
    compute_valence_data(operators, core_mask_list, valence_mask_list, method)

Compute the projected valence-only data for each operator using the given projection
`method` (`<:AbstractBasisProjection`). All operators must be hermitian and exactly one 
must be an `Overlap`. Returns a vector of wrapped `DenseRecipData` (one per input operator).
"""
function compute_valence_data(
    operators, core_mask_list, valence_mask_list, method::LaikovCore
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

    # Precompute core-core overlap inverse
    S₁₁ = unwrap_data(op_data(overlap))[core_mask_overlap, core_mask_overlap]
    S₁₂ = unwrap_data(op_data(overlap))[core_mask_overlap, valence_mask_overlap]
    S₁₁⁻¹ = inv(S₁₁)

    M_overlap = typeof(op_metadata(overlap))
    overlap_data₁₂ = wrap_data(M_overlap, S₁₂)
    overlap_data₁₁⁻¹ = wrap_data(M_overlap, S₁₁⁻¹)

    return [
        compute_valence_data(
            typeof(op_metadata(operator)),
            op_data(operator),
            overlap_data₁₁⁻¹,
            overlap_data₁₂,
            core_mask,
            valence_mask,
            method,
        )
        for (operator, core_mask, valence_mask) in
        zip(operators, core_mask_list, valence_mask_list)
    ]
end

# TODO: could use the same function if used abstract type of laikov and FC99V
# but would have to move Ŝ₁₂ outside
function compute_valence_data(
    ::Type{M},
    operator_data::DenseRecipData{T},
    overlap_data₁₁⁻¹::DenseRecipData{T},
    overlap_data₁₂::DenseRecipData{T},
    core_mask::BitVector,
    valence_mask::BitVector,
    ::LaikovCore,
) where {M<:AbstractMetadata,T}
    S₁₁⁻¹ = unwrap_data(overlap_data₁₁⁻¹)
    S₁₂ = unwrap_data(overlap_data₁₂)

    O₁₁ = unwrap_data(operator_data)[core_mask, core_mask]
    O₁₂ = unwrap_data(operator_data)[core_mask, valence_mask]
    O₂₂ = unwrap_data(operator_data)[valence_mask, valence_mask]

    # Ŝ₁₂ = S₁₁⁻¹ * S₁₂
    Ŝ₁₂ = LinearAlgebra.BLAS.hemm('L', 'U', S₁₁⁻¹, S₁₂)

    α = T(0.5)
    β = T(-1.0)
    LinearAlgebra.BLAS.hemm!('L', 'U', α, O₁₁, Ŝ₁₂, β, O₁₂)

    # A₂₂ = Ŝ₁₂' * (0.5 * O₁₁ * Ŝ₁₂ - O₁₂)
    A₂₂ = LinearAlgebra.BLAS.gemm('C', 'N', Ŝ₁₂, O₁₂)

    Ō₂₂ = O₂₂ + A₂₂ + A₂₂'

    return wrap_data(M, Ō₂₂)
end

### PARSING INTERFACE ###

get_basis_projection(::Val{:laikovcore}) = LaikovCore