struct LaikovCore <: AbstractBasisProjection end

function compute_valence_operator_data(operators, core_mask_list, valence_mask_list, method::LaikovCore)
    # Check that all operators are hermitian
    @argcheck all(operator -> get_sparsity(operator).hermitian, operators)

    # Requires exactly one overlap matrix for computing the projection
    # (if more than two than it's ambiguous which one to use)
    i_overlap = findall(operator -> get_kind(operator) isa Overlap, operators)
    @argcheck length(i_overlap) == 1

    overlap = operators[first(i_overlap)]
    core_mask_overlap = core_mask_list[first(i_overlap)]
    valence_mask_overlap = valence_mask_list[first(i_overlap)]

    # Precompute core-core overlap inverse
    S₁₁, S₁₂, _ = core_valence_partition(overlap, core_mask_overlap, valence_mask_overlap)
    S₁₁⁻¹ = inv(S₁₁)
    
    return [
        compute_valence_operator_data(operator, S₁₁⁻¹, S₁₂, core_mask, valence_mask, method)
        for (operator, core_mask, valence_mask) in zip(operators, core_mask_list, valence_mask_list)
    ]
end

function compute_valence_operator_data(operator::DenseOperator, S₁₁⁻¹::AbstractMatrix, S₁₂::AbstractMatrix,
    core_mask::BitVector, valence_mask::BitVector, ::LaikovCore)

    O₁₁, O₁₂, O₂₂ = core_valence_partition(operator, core_mask, valence_mask)

    Ŝ₁₂ = S₁₁⁻¹ * S₁₂
    A₂₂ = (0.5 * Ŝ₁₂' * O₁₁ - O₁₂') * Ŝ₁₂
    Ō₂₂ = O₂₂ + A₂₂ + A₂₂'

    return Ō₂₂
end

### PARSING INTERFACE ###

get_basis_projection(::Val{:laikovcore}) = LaikovCore