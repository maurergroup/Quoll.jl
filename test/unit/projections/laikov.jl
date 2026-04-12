#! format: off
using Quoll
using StaticArrays
using Dictionaries
using AtomsBase
using Unitful

using Test
using Main.TestFixtures

### high-level compute_valence_data dispatch ###

@testset "LaikovCore high-level dispatch" begin

    @testset "requires hermitian operators" begin
        overlap = make_canonical_dense_recip_operator(;
            hermitian=false, kind=Quoll.Overlap(; source=:ref),
            value=ComplexF64(1.0),
        )
        hamiltonian = make_canonical_dense_recip_operator(;
            hermitian=false, kind=Quoll.Hamiltonian(; source=:ref),
            value=ComplexF64(1.0),
        )
        operators = [hamiltonian, overlap]

        core_basis = make_test_core_basis()
        basisset = Quoll.op_basisset(overlap)
        core_masks = [
            Quoll.get_dense_subbasis_mask(basisset, core_basis; inverted=false)
            for _ in operators
        ]
        valence_masks = [
            Quoll.get_dense_subbasis_mask(basisset, core_basis; inverted=true)
            for _ in operators
        ]

        @test_throws ArgumentError Quoll.Projections.compute_valence_data(
            operators, core_masks, valence_masks, Quoll.Projections.LaikovCore()
        )
    end

    @testset "requires exactly one overlap (zero overlaps)" begin
        h1 = make_hermitian_dense_recip_operator(;
            kind=Quoll.Hamiltonian(; source=:ref),
        )
        h2 = make_hermitian_dense_recip_operator(;
            kind=Quoll.Hamiltonian(; source=:ref),
        )
        operators = [h1, h2]

        core_basis = make_test_core_basis()
        basisset = Quoll.op_basisset(h1)
        core_masks = [
            Quoll.get_dense_subbasis_mask(basisset, core_basis; inverted=false)
            for _ in operators
        ]
        valence_masks = [
            Quoll.get_dense_subbasis_mask(basisset, core_basis; inverted=true)
            for _ in operators
        ]

        @test_throws ArgumentError Quoll.Projections.compute_valence_data(
            operators, core_masks, valence_masks, Quoll.Projections.LaikovCore()
        )
    end

    @testset "requires exactly one overlap (two overlaps)" begin
        s1 = make_hermitian_dense_recip_operator(;
            kind=Quoll.Overlap(; source=:ref),
        )
        s2 = make_hermitian_dense_recip_operator(;
            kind=Quoll.Overlap(; source=:ref),
        )
        operators = [s1, s2]

        core_basis = make_test_core_basis()
        basisset = Quoll.op_basisset(s1)
        core_masks = [
            Quoll.get_dense_subbasis_mask(basisset, core_basis; inverted=false)
            for _ in operators
        ]
        valence_masks = [
            Quoll.get_dense_subbasis_mask(basisset, core_basis; inverted=true)
            for _ in operators
        ]

        @test_throws ArgumentError Quoll.Projections.compute_valence_data(
            operators, core_masks, valence_masks, Quoll.Projections.LaikovCore()
        )
    end

    @testset "nominal case" begin
        overlap = make_hermitian_dense_recip_operator(;
            kind=Quoll.Overlap(; source=:ref),
        )
        hamiltonian = make_hermitian_dense_recip_operator(;
            kind=Quoll.Hamiltonian(; source=:ref),
        )
        operators = [hamiltonian, overlap]

        core_basis = make_test_core_basis()
        basisset = Quoll.op_basisset(overlap)
        core_masks = [
            Quoll.get_dense_subbasis_mask(basisset, core_basis; inverted=false)
            for _ in operators
        ]
        valence_masks = [
            Quoll.get_dense_subbasis_mask(basisset, core_basis; inverted=true)
            for _ in operators
        ]

        results = Quoll.Projections.compute_valence_data(
            operators, core_masks, valence_masks, Quoll.Projections.LaikovCore()
        )

        @test length(results) == length(operators)

        n_valence = count(valence_masks[1])
        for result in results
            body = Quoll.unwrap_data(result)
            @test size(body) == (n_valence, n_valence)
            @test result isa Quoll.DenseRecipData
        end
    end

end

### per-operator compute_valence_data ###

@testset "LaikovCore per-operator compute_valence_data" begin

    @testset "output is DenseRecipData with correct dimensions" begin
        operator = make_hermitian_dense_recip_operator(;
            kind=Quoll.Hamiltonian(; source=:ref),
        )
        overlap = make_hermitian_dense_recip_operator(;
            kind=Quoll.Overlap(; source=:ref),
        )

        core_basis = make_test_core_basis()
        basisset = Quoll.op_basisset(operator)
        core_mask = Quoll.get_dense_subbasis_mask(basisset, core_basis; inverted=false)
        valence_mask = Quoll.get_dense_subbasis_mask(basisset, core_basis; inverted=true)

        # Prepare overlap data for per-operator call
        S₁₁ = Quoll.unwrap_data(Quoll.op_data(overlap))[core_mask, core_mask]
        S₁₂ = Quoll.unwrap_data(Quoll.op_data(overlap))[core_mask, valence_mask]
        S₁₁⁻¹ = inv(S₁₁)

        M = typeof(Quoll.op_metadata(overlap))
        overlap_data₁₂ = Quoll.wrap_data(M, S₁₂)
        overlap_data₁₁⁻¹ = Quoll.wrap_data(M, S₁₁⁻¹)

        result = Quoll.Projections.compute_valence_data(
            typeof(Quoll.op_metadata(operator)),
            Quoll.op_data(operator),
            overlap_data₁₁⁻¹,
            overlap_data₁₂,
            core_mask,
            valence_mask,
            Quoll.Projections.LaikovCore(),
        )

        @test result isa Quoll.DenseRecipData
        n_valence = count(valence_mask)
        @test size(Quoll.unwrap_data(result)) == (n_valence, n_valence)
    end

end
