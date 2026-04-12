#! format: off
using Quoll
using StaticArrays
using Dictionaries
using AtomsBase
using Unitful

using Test
using Main.TestFixtures

### fourier_transform_data! dispatch ###

@testset "fourier_transform_data! dispatch" begin

    @testset "Canonical KeyedOp → Recip Op (HasKeydata → NoKeydata)" begin
        # Input: Canonical KeyedOperator (HasKeydata, real-space) with nonzero data
        in_operator = make_canonical_block_real_keyed_operator(; value=1.0)

        # Output: Canonical DenseRecip Operator (NoKeydata, reciprocal-space)
        kpoint = SA[0.0, 0.0, 0.0]
        out_operator = make_canonical_dense_recip_operator(; kpoint=kpoint)

        # Compute phases for the Fourier transform
        images = Quoll.op_images(Quoll.op_sparsity(in_operator))
        phases_k = Quoll.precompute_phases([kpoint], images)[:, 1]

        # Dispatch: (NoKeydata, HasKeydata) → (Dₒᵤₜ, KDᵢₙ, Dᵢₙ)
        Quoll.fourier_transform_data!(out_operator, in_operator, phases_k)

        @test Quoll.op_data(out_operator) isa Quoll.DenseRecipData

        # With k=0 and phase=1, the FT sums over images. With all-ones input,
        # some output elements should be nonzero
        out_body = Quoll.unwrap_data(Quoll.op_data(out_operator))
        @test any(out_body .!= 0.0)
    end

    @testset "output is complex" begin
        in_operator = make_canonical_block_real_keyed_operator(; value=1.0)
        kpoint = SA[0.25, 0.0, 0.0]
        out_operator = make_canonical_dense_recip_operator(; kpoint=kpoint)

        images = Quoll.op_images(Quoll.op_sparsity(in_operator))
        phases_k = Quoll.precompute_phases([kpoint], images)[:, 1]

        Quoll.fourier_transform_data!(out_operator, in_operator, phases_k)

        out_body = Quoll.unwrap_data(Quoll.op_data(out_operator))
        @test eltype(out_body) <: Complex
    end

end

### inv_fourier_transform_data! dispatch ###

@testset "inv_fourier_transform_data! dispatch" begin

    @testset "Recip Op → Canonical KeyedOp (NoKeydata → HasKeydata)" begin
        # Input: Canonical DenseRecip Operator (NoKeydata) with nonzero data
        kpoint = SA[0.0, 0.0, 0.0]
        in_metadata = make_canonical_dense_recip_metadata(; kpoint=kpoint)
        in_operator = Quoll.build_operator(
            Quoll.Operator, in_metadata; value=Complex{Float64}(1.0)
        )

        # Output: Canonical KeyedOperator (HasKeydata, real-space)
        out_operator = make_canonical_block_real_keyed_operator()

        images = Quoll.op_images(Quoll.op_sparsity(out_operator))
        phases_k = Quoll.precompute_phases([kpoint], images)[:, 1]
        weight = 1.0

        # Dispatch: (HasKeydata, NoKeydata) → (KDₒᵤₜ, Dₒᵤₜ, Dᵢₙ)
        Quoll.inv_fourier_transform_data!(out_operator, in_operator, phases_k, weight)

        @test Quoll.op_keydata(out_operator) isa Quoll.CanonicalBlockRealKeyData

        # Data should have been filled by inverse FT
        keydata_body = Quoll.unwrap_data(Quoll.op_keydata(out_operator))
        has_nonzero = any(block -> any(block .!= 0.0), values(keydata_body))
        @test has_nonzero
    end

end

### fourier_transform high-level pipeline ###

@testset "fourier_transform pipeline" begin

    @testset "single kpoint" begin
        in_operator = make_canonical_block_real_keyed_operator(; value=1.0)
        kpoint = SA[0.0, 0.0, 0.0]

        out_operators = Quoll.fourier_transform(
            Quoll.Operator, Quoll.CanonicalDenseRecipMetadata,
            in_operator, [kpoint];
        )

        @test length(out_operators) == 1
        out_op = out_operators[1]

        @test out_op isa Quoll.Operator
        @test Quoll.op_source(out_op) isa Quoll.CanonicalSource
        @test Quoll.op_sparsity(out_op) isa Quoll.DenseRecipSparsity
        @test Quoll.op_data(out_op) isa Quoll.DenseRecipData
        @test Quoll.op_kpoint(out_op) ≈ kpoint

        # Output should be complex
        out_body = Quoll.unwrap_data(Quoll.op_data(out_op))
        @test eltype(out_body) <: Complex
        @test any(out_body .!= 0.0)
    end

    @testset "multiple kpoints" begin
        in_operator = make_canonical_block_real_keyed_operator(; value=1.0)
        kpoints = [SA[0.0, 0.0, 0.0], SA[0.25, 0.0, 0.0], SA[0.5, 0.5, 0.0]]

        out_operators = Quoll.fourier_transform(
            Quoll.Operator, Quoll.CanonicalDenseRecipMetadata,
            in_operator, kpoints;
        )

        @test length(out_operators) == 3

        for (out_op, kpt) in zip(out_operators, kpoints)
            @test out_op isa Quoll.Operator
            @test Quoll.op_kpoint(out_op) ≈ kpt
            @test Quoll.op_data(out_op) isa Quoll.DenseRecipData
        end
    end

    @testset "metadata preserved through fourier pipeline" begin
        in_operator = make_canonical_block_real_keyed_operator(; value=1.0)
        kpoint = SA[0.1, 0.2, 0.3]

        out_ops = Quoll.fourier_transform(
            Quoll.Operator, Quoll.CanonicalDenseRecipMetadata,
            in_operator, [kpoint];
        )

        out_metadata = Quoll.op_metadata(out_ops[1])
        @test out_metadata isa Quoll.RecipMetadata
        @test Quoll.op_source(out_metadata) isa Quoll.CanonicalSource

        # Basisset should be preserved (same source, no SH change)
        in_basisset = Quoll.op_basisset(in_operator)
        out_basisset = Quoll.op_basisset(out_ops[1])
        for z in keys(in_basisset.basis)
            @test in_basisset.basis[z] == out_basisset.basis[z]
        end
    end

    @testset "data dimensions match total basis size" begin
        in_operator = make_canonical_block_real_keyed_operator(; value=1.0)
        kpoint = SA[0.0, 0.0, 0.0]

        out_ops = Quoll.fourier_transform(
            Quoll.Operator, Quoll.CanonicalDenseRecipMetadata,
            in_operator, [kpoint];
        )

        basisset = Quoll.op_basisset(out_ops[1])
        nb = sum(Quoll.get_atom2nbasis(basisset))
        out_body = Quoll.unwrap_data(Quoll.op_data(out_ops[1]))
        @test size(out_body) == (nb, nb)
    end

end

### fourier roundtrip structure ###

@testset "fourier roundtrip structure" begin

    @testset "FT then inverse FT preserves structure" begin
        in_operator = make_canonical_block_real_keyed_operator(; value=1.0)
        kpoint = SA[0.0, 0.0, 0.0]

        # Forward FT: real-space KeyedOp → reciprocal Operator
        recip_ops = Quoll.fourier_transform(
            Quoll.Operator, Quoll.CanonicalDenseRecipMetadata,
            in_operator, [kpoint];
        )
        recip_op = recip_ops[1]

        # Inverse FT: reciprocal Operator → real-space KeyedOp
        out_metadata = make_canonical_block_real_metadata()
        out_operator = Quoll.build_operator(Quoll.KeyedOperator, out_metadata)

        images = Quoll.op_images(Quoll.op_sparsity(out_operator))
        phases_k = Quoll.precompute_phases([kpoint], images)[:, 1]
        weight = 1.0

        Quoll.inv_fourier_transform_data!(out_operator, recip_op, phases_k, weight)

        # Structure should be preserved
        @test out_operator isa Quoll.KeyedOperator
        @test Quoll.op_source(out_operator) isa Quoll.CanonicalSource
        @test Quoll.op_data(out_operator) isa Quoll.CanonicalBlockRealData
        @test Quoll.op_keydata(out_operator) isa Quoll.CanonicalBlockRealKeyData

        # Data should be nonzero after roundtrip
        keydata_body = Quoll.unwrap_data(Quoll.op_keydata(out_operator))
        has_nonzero = any(block -> any(block .!= 0.0), values(keydata_body))
        @test has_nonzero
    end

end
