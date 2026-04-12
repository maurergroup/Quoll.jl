#! format: off
using Quoll
using StaticArrays
using Dictionaries
using AtomsBase
using Unitful

using Test
using Main.TestFixtures

@testset "convert_sparsity (metadata-level)" begin

    @testset "same type, preserve hermicity" begin
        metadata = make_canonical_block_real_metadata(; hermitian=false)
        out = Quoll.convert_sparsity(Quoll.BlockRealSparsity, metadata; hermitian=false)

        @test out isa Quoll.BlockRealSparsity
        @test out.hermitian == false
        @test sort(collect(keys(out.ij2images))) == sort(collect(keys(Quoll.op_sparsity(metadata).ij2images)))
    end

    @testset "nonhermitian to hermitian" begin
        metadata = make_canonical_block_real_metadata(; hermitian=false)
        out = Quoll.convert_sparsity(Quoll.BlockRealSparsity, metadata; hermitian=true)

        @test out isa Quoll.BlockRealSparsity
        @test out.hermitian == true
        # Hermitian keeps only upper triangle: i <= j
        for (i, j) in keys(out.ij2images)
            @test i <= j
        end
    end

    @testset "hermitian to nonhermitian" begin
        metadata = make_canonical_block_real_metadata(; hermitian=true)
        out = Quoll.convert_sparsity(Quoll.BlockRealSparsity, metadata; hermitian=false)

        @test out isa Quoll.BlockRealSparsity
        @test out.hermitian == false
        # Nonhermitian should have mirror pairs for off-diagonal
        has_12 = (1, 2) ∈ keys(out.ij2images)
        has_21 = (2, 1) ∈ keys(out.ij2images)
        @test has_12
        @test has_21
    end

    @testset "hermitian defaults to metadata" begin
        metadata_nh = make_canonical_block_real_metadata(; hermitian=false)
        out_nh = Quoll.convert_sparsity(Quoll.BlockRealSparsity, metadata_nh; hermitian=nothing)
        @test out_nh.hermitian == false

        metadata_h = make_canonical_block_real_metadata(; hermitian=true)
        out_h = Quoll.convert_sparsity(Quoll.BlockRealSparsity, metadata_h; hermitian=nothing)
        @test out_h.hermitian == true
    end

    @testset "with radii" begin
        metadata = make_canonical_block_real_metadata(; hermitian=false)
        # Use a small radius that should produce a different sparsity than input
        radii = Dict(
            ChemicalSpecies(:H)  => 0.5u"Å",
            ChemicalSpecies(:Li) => 0.5u"Å",
        )
        out = Quoll.convert_sparsity(
            Quoll.BlockRealSparsity, metadata; radii=radii, hermitian=false
        )

        @test out isa Quoll.BlockRealSparsity
        @test out.hermitian == false
        # With small radii and atoms far apart, only onsite terms should be present
        for (i, j) in keys(out.ij2images)
            @test i == j
        end
    end

    @testset "BlockReal to DenseRecip" begin
        metadata = make_canonical_block_real_metadata(; hermitian=false)
        out = Quoll.convert_sparsity(Quoll.DenseRecipSparsity, metadata; hermitian=false)

        @test out isa Quoll.DenseRecipSparsity
        @test out.hermitian == false
    end

    @testset "BlockReal to DenseRecip hermitian" begin
        metadata = make_canonical_block_real_metadata(; hermitian=false)
        out = Quoll.convert_sparsity(Quoll.DenseRecipSparsity, metadata; hermitian=true)

        @test out isa Quoll.DenseRecipSparsity
        @test out.hermitian == true
    end

end
