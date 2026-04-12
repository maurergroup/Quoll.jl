#! format: off
using Quoll
using StaticArrays
using Dictionaries
using AtomsBase
using Unitful

using Test
using Main.TestFixtures

@testset "core_valence_partition" begin

    operator = make_hermitian_dense_recip_operator()

    core_basis = make_test_core_basis()
    basisset = Quoll.op_basisset(operator)
    core_mask = Quoll.get_dense_subbasis_mask(basisset, core_basis; inverted=false)
    valence_mask = Quoll.get_dense_subbasis_mask(basisset, core_basis; inverted=true)

    O₁₁, O₁₂, O₂₂ = Quoll.Projections.core_valence_partition(
        Quoll.op_data(operator), core_mask, valence_mask
    )

    n_core = count(core_mask)
    n_valence = count(valence_mask)

    @test size(O₁₁) == (n_core, n_core)
    @test size(O₁₂) == (n_core, n_valence)
    @test size(O₂₂) == (n_valence, n_valence)

    @test O₁₁ isa SubArray
    @test O₁₂ isa SubArray
    @test O₂₂ isa SubArray

end
