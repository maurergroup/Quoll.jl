#! format: off
using Quoll
using StaticArrays
using Dictionaries
using AtomsBase
using Unitful

using Test
using Main.TestFixtures

### convert_data! dispatch ###

@testset "convert_data! dispatch" begin

    @testset "DeepH Op → Canonical KeyedOp (NoKeydata → HasKeydata)" begin
        # Input: DeepH Operator (NoKeydata) with nonzero data
        in_metadata = make_deeph_block_real_metadata(; hermitian=false)
        in_operator = Quoll.build_operator(Quoll.Operator, in_metadata; value=1.0)

        # Output: Canonical KeyedOperator (HasKeydata) with matching metadata
        out_metadata = Quoll.convert_metadata(Quoll.CanonicalBlockRealMetadata, in_metadata)
        out_operator = Quoll.build_operator(Quoll.KeyedOperator, out_metadata)

        # Dispatch should route through (HasKeydata, NoKeydata) → (KDₒᵤₜ, Dₒᵤₜ, Dᵢₙ)
        Quoll.convert_data!(out_operator, in_operator)

        @test out_operator isa Quoll.KeyedOperator
        @test Quoll.trait(Quoll.KeyedTrait, typeof(out_operator)) == Quoll.HasKeydata()
        @test Quoll.op_data(out_operator) isa Quoll.CanonicalBlockRealData
        @test Quoll.op_keydata(out_operator) isa Quoll.CanonicalBlockRealKeyData

        # Data should have been transferred (not all zeros)
        keydata_body = Quoll.unwrap_data(Quoll.op_keydata(out_operator))
        has_nonzero = any(block -> any(block .!= 0.0), values(keydata_body))
        @test has_nonzero
    end

    @testset "Canonical KeyedOp → DeepH Op (HasKeydata → NoKeydata)" begin
        # Input: Canonical KeyedOperator (HasKeydata) with nonzero data
        in_metadata = make_canonical_block_real_metadata(; hermitian=false)
        in_operator = Quoll.build_operator(Quoll.KeyedOperator, in_metadata; value=1.0)

        # Output: DeepH Operator (NoKeydata) with matching metadata
        out_metadata = Quoll.convert_metadata(Quoll.DeepHBlockRealMetadata, in_metadata)
        out_operator = Quoll.build_operator(Quoll.Operator, out_metadata)

        # Dispatch should route through (NoKeydata, HasKeydata) → (Dₒᵤₜ, KDᵢₙ, Dᵢₙ)
        Quoll.convert_data!(out_operator, in_operator)

        @test out_operator isa Quoll.Operator
        @test Quoll.trait(Quoll.KeyedTrait, typeof(out_operator)) == Quoll.NoKeydata()
        @test Quoll.op_data(out_operator) isa Quoll.DeepHBlockRealData

        # Data should have been transferred
        data_body = Quoll.unwrap_data(Quoll.op_data(out_operator))
        has_nonzero = any(block -> any(block .!= 0.0), values(data_body))
        @test has_nonzero
    end

    @testset "output data types match metadata" begin
        # DeepH → Canonical: verify data type aliases are correctly resolved by dispatch
        in_metadata = make_deeph_block_real_metadata(; hermitian=false)
        in_op = Quoll.build_operator(Quoll.Operator, in_metadata; value=1.0)

        out_metadata = Quoll.convert_metadata(Quoll.CanonicalBlockRealMetadata, in_metadata)
        out_op = Quoll.build_operator(Quoll.KeyedOperator, out_metadata)
        Quoll.convert_data!(out_op, in_op)

        Mₒᵤₜ = typeof(Quoll.op_metadata(out_op))
        @test typeof(Quoll.op_data(out_op)) <: Quoll.op_data_type(Mₒᵤₜ)
        @test typeof(Quoll.op_keydata(out_op)) <: Quoll.op_keydata_type(Mₒᵤₜ)
    end

end

### convert_operator end-to-end ###

@testset "convert_operator" begin

    @testset "Canonical KeyedOp → DeepH Op" begin
        in_operator = make_canonical_block_real_keyed_operator(; value=1.0)

        out_operator = Quoll.convert_operator(
            Quoll.Operator, Quoll.DeepHBlockRealMetadata, in_operator
        )

        # Verify full pipeline result
        @test out_operator isa Quoll.Operator
        @test Quoll.op_source(out_operator) isa Quoll.DeepHSource
        @test Quoll.op_sparsity(out_operator) isa Quoll.BlockRealSparsity
        @test Quoll.op_shconv(out_operator) == Quoll.default_shconv(Quoll.DeepHSource())
        @test Quoll.op_data(out_operator) isa Quoll.DeepHBlockRealData

        # Data should have been converted
        data_body = Quoll.unwrap_data(Quoll.op_data(out_operator))
        has_nonzero = any(block -> any(block .!= 0.0), values(data_body))
        @test has_nonzero
    end

    @testset "DeepH Op → Canonical KeyedOp" begin
        in_operator = make_deeph_block_real_operator(; value=1.0)

        out_operator = Quoll.convert_operator(
            Quoll.KeyedOperator, Quoll.CanonicalBlockRealMetadata, in_operator
        )

        @test out_operator isa Quoll.KeyedOperator
        @test Quoll.op_source(out_operator) isa Quoll.CanonicalSource
        @test Quoll.op_sparsity(out_operator) isa Quoll.BlockRealSparsity
        @test Quoll.op_shconv(out_operator) == Quoll.default_shconv(Quoll.CanonicalSource())
        @test Quoll.op_data(out_operator) isa Quoll.CanonicalBlockRealData
        @test Quoll.op_keydata(out_operator) isa Quoll.CanonicalBlockRealKeyData

        keydata_body = Quoll.unwrap_data(Quoll.op_keydata(out_operator))
        has_nonzero = any(block -> any(block .!= 0.0), values(keydata_body))
        @test has_nonzero
    end

    @testset "metadata is correctly converted in pipeline" begin
        in_operator = make_deeph_block_real_operator(; value=1.0)

        out_operator = Quoll.convert_operator(
            Quoll.KeyedOperator, Quoll.CanonicalBlockRealMetadata, in_operator
        )

        # The output metadata should be a proper Canonical metadata
        out_metadata = Quoll.op_metadata(out_operator)
        @test out_metadata isa Quoll.RealMetadata
        @test Quoll.op_source(out_metadata) isa Quoll.CanonicalSource

        # Basisset should be reordered to Canonical SH convention
        H = ChemicalSpecies(:H)
        p_orbs = filter(b -> b.l == 1, Quoll.op_basisset(out_metadata).basis[H])
        @test [b.m for b in p_orbs] == [-1, 0, 1]
    end

    @testset "hermicity preserved through pipeline" begin
        in_operator = make_canonical_block_real_keyed_operator(; hermitian=false, value=1.0)
        out_operator = Quoll.convert_operator(
            Quoll.Operator, Quoll.DeepHBlockRealMetadata, in_operator
        )

        @test Quoll.op_hermicity(out_operator) == false
    end

    @testset "float type preserved through pipeline" begin
        in_metadata = make_canonical_block_real_metadata(; hermitian=false)
        in_operator = Quoll.build_operator(
            Quoll.KeyedOperator, in_metadata; value=Float32(1.0)
        )

        out_operator = Quoll.convert_operator(
            Quoll.Operator, Quoll.DeepHBlockRealMetadata, in_operator
        )

        out_data = Quoll.op_data(out_operator)
        first_block = first(values(Quoll.unwrap_data(out_data)))
        @test eltype(first_block) == Float32
    end

    @testset "explicit hermitian kwarg" begin
        in_operator = make_deeph_block_real_operator(; hermitian=false, value=1.0)

        out_operator = Quoll.convert_operator(
            Quoll.KeyedOperator, Quoll.CanonicalBlockRealMetadata, in_operator;
            hermitian=true,
        )

        @test Quoll.op_hermicity(out_operator) == true
    end

end

### roundtrip structure preservation ###

@testset "convert_data! roundtrip structure" begin

    @testset "Canonical → DeepH → Canonical preserves structure" begin
        orig = make_canonical_block_real_keyed_operator(; value=2.0)

        # Canonical → DeepH
        mid_metadata = Quoll.convert_metadata(
            Quoll.DeepHBlockRealMetadata, Quoll.op_metadata(orig)
        )
        mid = Quoll.build_operator(Quoll.Operator, mid_metadata)
        Quoll.convert_data!(mid, orig)

        # DeepH → Canonical
        final_metadata = Quoll.convert_metadata(
            Quoll.CanonicalBlockRealMetadata, Quoll.op_metadata(mid)
        )
        final = Quoll.build_operator(Quoll.KeyedOperator, final_metadata)
        Quoll.convert_data!(final, mid)

        # Structure should be preserved
        @test final isa Quoll.KeyedOperator
        @test Quoll.op_source(final) isa Quoll.CanonicalSource
        @test Quoll.op_data(final) isa Quoll.CanonicalBlockRealData
        @test Quoll.op_keydata(final) isa Quoll.CanonicalBlockRealKeyData

        # Data should survive roundtrip (nonzero)
        keydata_body = Quoll.unwrap_data(Quoll.op_keydata(final))
        has_nonzero = any(block -> any(block .!= 0.0), values(keydata_body))
        @test has_nonzero
    end

end
