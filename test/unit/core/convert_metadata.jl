#! format: off
using Quoll
using StaticArrays
using Dictionaries
using AtomsBase
using Unitful

using Test
using Main.TestFixtures

### convert_metadata_basic ###

@testset "convert_metadata_basic" begin

    @testset "same format identity" begin
        in_metadata = make_canonical_block_real_metadata(; hermitian=false)
        out = Quoll.convert_metadata_basic(
            Quoll.CanonicalBlockRealMetadata, in_metadata;
        )

        @test out isa Quoll.BasicMetadataContainer
        @test Quoll.op_source(out) isa Quoll.CanonicalSource
        @test Quoll.op_sparsity(out) isa Quoll.BlockRealSparsity
        @test Quoll.op_shconv(out) == Quoll.default_shconv(Quoll.CanonicalSource())

        # Basisset ordering should be preserved (Canonical → Canonical is identity SH change)
        in_basis = Quoll.op_basisset(in_metadata).basis
        out_basis = Quoll.op_basisset(out).basis
        for z in keys(in_basis)
            @test in_basis[z] == out_basis[z]
        end
    end

    @testset "cross-format SH reorder (DeepH to Canonical)" begin
        in_metadata = make_deeph_block_real_metadata(; hermitian=false)
        out = Quoll.convert_metadata_basic(
            Quoll.CanonicalBlockRealMetadata, in_metadata;
        )

        @test Quoll.op_source(out) isa Quoll.CanonicalSource
        @test Quoll.op_shconv(out) == Quoll.default_shconv(Quoll.CanonicalSource())

        # DeepH p ordering is [p_1, p_{-1}, p_0]. Canonical is [p_{-1}, p_0, p_1].
        # After conversion, basis should be back in wiki (canonical) order.
        H = ChemicalSpecies(:H)
        out_basis_H = Quoll.op_basisset(out).basis[H]
        p_orbs = filter(b -> b.l == 1, out_basis_H)
        @test [b.m for b in p_orbs] == [-1, 0, 1]
    end

    @testset "custom out_shconv" begin
        in_metadata = make_canonical_block_real_metadata(; hermitian=false)
        custom_shconv = Quoll.default_shconv(Quoll.DeepHSource())

        out = Quoll.convert_metadata_basic(
            Quoll.CanonicalBlockRealMetadata, in_metadata;
            out_shconv=custom_shconv,
        )

        @test Quoll.op_shconv(out) == custom_shconv
    end

    @testset "with subbasis (keep only p orbitals)" begin
        in_metadata = make_canonical_block_real_metadata(; hermitian=false)

        H = ChemicalSpecies(:H)
        Li = ChemicalSpecies(:Li)
        subbasis = [
            Quoll.BasisMetadata(H,  2, 1, -1, nothing),
            Quoll.BasisMetadata(H,  2, 1,  0, nothing),
            Quoll.BasisMetadata(H,  2, 1,  1, nothing),
            Quoll.BasisMetadata(Li, 2, 1, -1, nothing),
            Quoll.BasisMetadata(Li, 2, 1,  0, nothing),
            Quoll.BasisMetadata(Li, 2, 1,  1, nothing),
        ]

        out = Quoll.convert_metadata_basic(
            Quoll.CanonicalBlockRealMetadata, in_metadata;
            subbasis=subbasis,
        )

        # Each species should have 3 basis functions (p only)
        out_basis = Quoll.op_basisset(out).basis
        @test length(out_basis[H]) == 3
        @test length(out_basis[Li]) == 3
        @test all(b -> b.l == 1, out_basis[H])
        @test all(b -> b.l == 1, out_basis[Li])
    end

    @testset "with subbasis inverted (remove p orbitals)" begin
        in_metadata = make_canonical_block_real_metadata(; hermitian=false)

        H = ChemicalSpecies(:H)
        Li = ChemicalSpecies(:Li)
        subbasis = [
            Quoll.BasisMetadata(H,  2, 1, -1, nothing),
            Quoll.BasisMetadata(H,  2, 1,  0, nothing),
            Quoll.BasisMetadata(H,  2, 1,  1, nothing),
            Quoll.BasisMetadata(Li, 2, 1, -1, nothing),
            Quoll.BasisMetadata(Li, 2, 1,  0, nothing),
            Quoll.BasisMetadata(Li, 2, 1,  1, nothing),
        ]

        out = Quoll.convert_metadata_basic(
            Quoll.CanonicalBlockRealMetadata, in_metadata;
            subbasis=subbasis, inverted=true,
        )

        # Each species should have 1 basis function (s only)
        out_basis = Quoll.op_basisset(out).basis
        @test length(out_basis[H]) == 1
        @test length(out_basis[Li]) == 1
        @test all(b -> b.l == 0, out_basis[H])
        @test all(b -> b.l == 0, out_basis[Li])
    end

    @testset "hermitian kwarg overrides default" begin
        in_metadata = make_canonical_block_real_metadata(; hermitian=false)
        out = Quoll.convert_metadata_basic(
            Quoll.CanonicalBlockRealMetadata, in_metadata;
            hermitian=true,
        )

        @test Quoll.op_hermicity(out) == true
    end

end

### convert_metadata_extra ###

@testset "convert_metadata_extra" begin

    @testset "RealMetadata has no extras" begin
        in_metadata = make_canonical_block_real_metadata()
        basic = Quoll.op_basic_metadata(in_metadata)

        extras = Quoll.convert_metadata_extra(
            Quoll.CanonicalBlockNoSpinRealMetadata, in_metadata, basic;
        )
        @test isempty(collect(extras))
    end

    @testset "RecipMetadata extracts kpoint from input" begin
        kpoint = SA[0.25, 0.5, 0.0]
        in_metadata = make_canonical_dense_recip_metadata(; kpoint=kpoint)
        basic = Quoll.op_basic_metadata(in_metadata)

        extras = Quoll.convert_metadata_extra(
            Quoll.CanonicalDenseRecipMetadata, in_metadata, basic;
        )
        extras_collected = collect(extras)
        @test length(extras_collected) == 1
        @test extras_collected[1] ≈ kpoint
    end

    @testset "Real to Recip without kpoint kwarg throws" begin
        in_metadata = make_canonical_block_real_metadata()
        basic = Quoll.op_basic_metadata(in_metadata)

        @test_throws ArgumentError Quoll.convert_metadata_extra(
            Quoll.CanonicalDenseRecipMetadata, in_metadata, basic;
        )
    end

    @testset "Real to Recip with kpoint kwarg" begin
        in_metadata = make_canonical_block_real_metadata()
        basic = Quoll.op_basic_metadata(in_metadata)
        kpoint = SA[0.1, 0.2, 0.3]

        extras = Quoll.convert_metadata_extra(
            Quoll.CanonicalDenseRecipMetadata, in_metadata, basic;
            extra_kwargs=(; kpoint=kpoint),
        )
        extras_collected = collect(extras)
        @test length(extras_collected) == 1
        @test extras_collected[1] ≈ kpoint
    end

    @testset "SpinRealMetadata passes spins through" begin
        in_metadata = make_spin_real_metadata(; spin=:up)
        basic = Quoll.op_basic_metadata(in_metadata)

        extras = Quoll.convert_metadata_extra(
            Quoll.CanonicalBlockSpinRealMetadata, in_metadata, basic;
        )
        extras_collected = collect(extras)
        @test length(extras_collected) == 1
        @test extras_collected[1] isa Quoll.SpinsMetadata
    end

    @testset "RealMetadata to SpinReal union returns missing for spin" begin
        in_metadata = make_canonical_block_real_metadata()
        basic = Quoll.op_basic_metadata(in_metadata)

        # Union type has [SpinTrait] in extrafield_traittypes, but NoSpin input
        # with RealMetadata input triggers the method returning missing
        extras = Quoll.convert_metadata_extra(
            Quoll.CanonicalBlockRealMetadata, in_metadata, basic;
        )
        extras_collected = collect(extras)
        @test isempty(extras_collected)
    end

end

### convert_metadata_final ###

@testset "convert_metadata_final" begin

    @testset "RealMetadata from union, no extras" begin
        metadata = make_canonical_block_real_metadata()
        basic = Quoll.op_basic_metadata(metadata)

        out = Quoll.convert_metadata_final(Quoll.CanonicalBlockRealMetadata, basic)
        @test out isa Quoll.RealMetadata
    end

    @testset "SpinRealMetadata from union, with spins" begin
        metadata = make_spin_real_metadata()
        basic = Quoll.op_basic_metadata(metadata)
        spins = Quoll.op_spins(metadata)

        out = Quoll.convert_metadata_final(Quoll.CanonicalBlockRealMetadata, basic, spins)
        @test out isa Quoll.SpinRealMetadata
        @test Quoll.op_spins(out) === spins
    end

    @testset "RecipMetadata from union, with kpoint" begin
        kpoint = SA[0.1, 0.2, 0.3]
        metadata = make_canonical_dense_recip_metadata(; kpoint=kpoint)
        basic = Quoll.op_basic_metadata(metadata)

        out = Quoll.convert_metadata_final(Quoll.CanonicalDenseRecipMetadata, basic, kpoint)
        @test out isa Quoll.RecipMetadata
        @test Quoll.op_kpoint(out) ≈ kpoint
    end

    @testset "SpinRecipMetadata from union, with kpoint and spins" begin
        kpoint = SA[0.1, 0.2, 0.3]
        metadata = make_spin_recip_metadata(; kpoint=kpoint)
        basic = Quoll.op_basic_metadata(metadata)
        spins = Quoll.op_spins(metadata)

        # Alphabetical trait sort: SpaceTrait < SpinTrait → (kpoint, spins)
        out = Quoll.convert_metadata_final(
            Quoll.CanonicalDenseRecipMetadata, basic, kpoint, spins
        )
        @test out isa Quoll.SpinRecipMetadata
        @test Quoll.op_kpoint(out) ≈ kpoint
        @test Quoll.op_spins(out) === spins
    end

end

### convert_metadata end-to-end ###

@testset "convert_metadata end-to-end" begin

    @testset "Canonical Real to Canonical Real (identity)" begin
        in_metadata = make_canonical_block_real_metadata(; hermitian=false)
        out = Quoll.convert_metadata(Quoll.CanonicalBlockRealMetadata, in_metadata)

        @test out isa Quoll.RealMetadata
        @test Quoll.op_source(out) isa Quoll.CanonicalSource
        @test Quoll.op_sparsity(out) isa Quoll.BlockRealSparsity
        @test Quoll.op_hermicity(out) == false
    end

    @testset "DeepH Real to Canonical Real" begin
        in_metadata = make_deeph_block_real_metadata(; hermitian=false)
        out = Quoll.convert_metadata(Quoll.CanonicalBlockRealMetadata, in_metadata)

        @test out isa Quoll.RealMetadata
        @test Quoll.op_source(out) isa Quoll.CanonicalSource
        @test Quoll.op_shconv(out) == Quoll.default_shconv(Quoll.CanonicalSource())

        # Verify basis is now in Canonical (wiki) ordering
        H = ChemicalSpecies(:H)
        p_orbs = filter(b -> b.l == 1, Quoll.op_basisset(out).basis[H])
        @test [b.m for b in p_orbs] == [-1, 0, 1]
    end

    @testset "Canonical Real to Canonical Recip (with kpoint)" begin
        in_metadata = make_canonical_block_real_metadata()
        kpoint = SA[0.25, 0.0, 0.0]

        out = Quoll.convert_metadata(
            Quoll.CanonicalDenseRecipMetadata, in_metadata;
            extra_kwargs=(; kpoint=kpoint),
        )

        @test out isa Quoll.RecipMetadata
        @test Quoll.op_source(out) isa Quoll.CanonicalSource
        @test Quoll.op_sparsity(out) isa Quoll.DenseRecipSparsity
        @test Quoll.op_kpoint(out) ≈ kpoint
    end

    @testset "SpinReal to SpinReal (spin preserved)" begin
        in_metadata = make_spin_real_metadata(; spin=:up)
        out = Quoll.convert_metadata(Quoll.CanonicalBlockRealMetadata, in_metadata)

        @test out isa Quoll.SpinRealMetadata
        @test Quoll.op_spins(out) isa Quoll.SpinsMetadata
    end

    @testset "SpinReal to SpinRecip (with kpoint)" begin
        in_metadata = make_spin_real_metadata(; spin=:up)
        kpoint = SA[0.5, 0.0, 0.0]

        out = Quoll.convert_metadata(
            Quoll.CanonicalDenseRecipMetadata, in_metadata;
            extra_kwargs=(; kpoint=kpoint),
        )

        @test out isa Quoll.SpinRecipMetadata
        @test Quoll.op_kpoint(out) ≈ kpoint
        @test Quoll.op_spins(out) isa Quoll.SpinsMetadata
    end

end
