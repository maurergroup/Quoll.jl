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
    S₁₁ = get_data(overlap)[core_mask_overlap, core_mask_overlap]
    S₁₂ = get_data(overlap)[core_mask_overlap, valence_mask_overlap]
    S₁₁⁻¹ = inv(S₁₁)
    
    return [
        compute_valence_operator_data(operator, S₁₁⁻¹, S₁₂, core_mask, valence_mask, method, get_float(operator))
        for (operator, core_mask, valence_mask) in zip(operators, core_mask_list, valence_mask_list)
    ]
end

# TODO: could use the same function if used abstract type of laikov and FC99V
# but would have to move Ŝ₁₂ outside
function compute_valence_operator_data(operator::DenseOperator, S₁₁⁻¹::AbstractMatrix, S₁₂::AbstractMatrix,
    core_mask::BitVector, valence_mask::BitVector, ::LaikovCore, ::Type{T}) where T<:Complex

    O₁₁ = get_data(operator)[core_mask, core_mask]
    O₁₂ = get_data(operator)[core_mask, valence_mask]
    O₂₂ = get_data(operator)[valence_mask, valence_mask]

    # Ŝ₁₂ = S₁₁⁻¹ * S₁₂
    Ŝ₁₂ = LinearAlgebra.BLAS.hemm('L', 'U', S₁₁⁻¹, S₁₂)

    α = T(0.5)
    β = T(-1.0)
    LinearAlgebra.BLAS.hemm!('L', 'U', α, O₁₁, Ŝ₁₂, β, O₁₂)

    # A₂₂ = Ŝ₁₂' * (0.5 * O₁₁ * Ŝ₁₂ - O₁₂)
    A₂₂ = LinearAlgebra.BLAS.gemm('C', 'N', Ŝ₁₂, O₁₂)

    Ō₂₂ = O₂₂ + A₂₂ + A₂₂'

    return Ō₂₂
end

### PARSING INTERFACE ###

get_basis_projection(::Val{:laikovcore}) = LaikovCore