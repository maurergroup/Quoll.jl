#! format: off
using Quoll
using StaticArrays
using Dictionaries
using AtomsBase
using Unitful

using Test
using Main.TestFixtures

@testset "build_operator" begin

    @testset "Operator, CanonicalBlockReal, default value" begin
        metadata = make_canonical_block_real_metadata(; hermitian=false)
        op = Quoll.build_operator(Quoll.Operator, metadata)

        @test op isa Quoll.Operator
        @test Quoll.op_metadata(op) === metadata
        @test Quoll.trait(Quoll.KeyedTrait, typeof(op)) == Quoll.NoKeydata()

        data = Quoll.op_data(op)
        @test data isa Quoll.CanonicalBlockRealData
        body = Quoll.unwrap_data(data)
        for arr in values(body)
            @test all(arr .== 0.0)
        end
    end

    @testset "Operator, CanonicalBlockReal, custom value" begin
        metadata = make_canonical_block_real_metadata(; hermitian=false)
        op = Quoll.build_operator(Quoll.Operator, metadata; value=1.5)

        body = Quoll.unwrap_data(Quoll.op_data(op))
        for arr in values(body)
            @test all(arr .== 1.5)
        end
    end

    @testset "Operator, CanonicalBlockReal, Float32 type" begin
        metadata = make_canonical_block_real_metadata(; hermitian=false)
        op = Quoll.build_operator(Quoll.Operator, metadata; value=Float32(0.0))

        body = Quoll.unwrap_data(Quoll.op_data(op))
        for arr in values(body)
            @test eltype(arr) == Float32
        end
    end

    @testset "Operator, CanonicalBlockReal, uninitialised" begin
        metadata = make_canonical_block_real_metadata(; hermitian=false)
        op = Quoll.build_operator(
            Quoll.Operator, metadata; type=Float64, initialised=false
        )

        # Just verify data is allocated with correct dimensions, not values
        body = Quoll.unwrap_data(Quoll.op_data(op))
        basisset = Quoll.op_basisset(metadata)
        species2nbasis = Quoll.get_species2nbasis(basisset)

        for ((z1, z2), arr) in pairs(body)
            @test size(arr, 1) == species2nbasis[z1]
            @test size(arr, 2) == species2nbasis[z2]
        end
    end

    @testset "Operator, DeepHBlockReal" begin
        metadata = make_deeph_block_real_metadata(; hermitian=false)
        op = Quoll.build_operator(Quoll.Operator, metadata)

        @test op isa Quoll.Operator
        data = Quoll.op_data(op)
        @test data isa Quoll.DeepHBlockRealData

        body = Quoll.unwrap_data(data)
        sparsity = Quoll.op_sparsity(metadata)
        ref_keys = Set(Quoll.get_deeph_keys(sparsity))
        @test Set(collect(keys(body))) == ref_keys

        for arr in values(body)
            @test all(arr .== 0.0)
        end
    end

    @testset "KeyedOperator, CanonicalBlockReal" begin
        metadata = make_canonical_block_real_metadata(; hermitian=false)
        op = Quoll.build_operator(Quoll.KeyedOperator, metadata)

        @test op isa Quoll.KeyedOperator
        @test Quoll.trait(Quoll.KeyedTrait, typeof(op)) == Quoll.HasKeydata()

        data = Quoll.op_data(op)
        keydata = Quoll.op_keydata(op)
        @test data isa Quoll.CanonicalBlockRealData
        @test keydata isa Quoll.CanonicalBlockRealKeyData

        # Keydata should have atom-pair keys matching sparsity
        sparsity = Quoll.op_sparsity(metadata)
        @test Set(collect(keys(Quoll.unwrap_data(keydata)))) == Set(collect(keys(sparsity.ij2images)))
    end

    @testset "Operator, CanonicalDenseRecip" begin
        metadata = make_canonical_dense_recip_metadata(; hermitian=false)
        op = Quoll.build_operator(Quoll.Operator, metadata)

        @test op isa Quoll.Operator
        data = Quoll.op_data(op)
        @test data isa Quoll.DenseRecipData

        body = Quoll.unwrap_data(data)
        basisset = Quoll.op_basisset(metadata)
        nb = sum(Quoll.get_atom2nbasis(basisset))
        @test size(body) == (nb, nb)
        @test all(body .== 0.0)
    end

    @testset "Operator with subbasis" begin
        metadata = make_canonical_block_real_metadata(; hermitian=false)

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

        op = Quoll.build_operator(Quoll.Operator, metadata; subbasis=subbasis)

        # Metadata should have reduced basisset (3 per species instead of 4)
        out_basisset = Quoll.op_basisset(op)
        @test length(out_basisset.basis[H]) == 3
        @test length(out_basisset.basis[Li]) == 3

        # Data dimensions should match reduced basis
        body = Quoll.unwrap_data(Quoll.op_data(op))
        species2nbasis = Quoll.get_species2nbasis(out_basisset)
        for ((z1, z2), arr) in pairs(body)
            @test size(arr, 1) == species2nbasis[z1]
            @test size(arr, 2) == species2nbasis[z2]
        end
    end

    @testset "Operator with subbasis inverted" begin
        metadata = make_canonical_block_real_metadata(; hermitian=false)

        H = ChemicalSpecies(:H)
        Li = ChemicalSpecies(:Li)
        # Remove p orbitals → keep only s
        subbasis = [
            Quoll.BasisMetadata(H,  2, 1, -1, nothing),
            Quoll.BasisMetadata(H,  2, 1,  0, nothing),
            Quoll.BasisMetadata(H,  2, 1,  1, nothing),
            Quoll.BasisMetadata(Li, 2, 1, -1, nothing),
            Quoll.BasisMetadata(Li, 2, 1,  0, nothing),
            Quoll.BasisMetadata(Li, 2, 1,  1, nothing),
        ]

        op = Quoll.build_operator(
            Quoll.Operator, metadata; subbasis=subbasis, inverted=true
        )

        out_basisset = Quoll.op_basisset(op)
        @test length(out_basisset.basis[H]) == 1
        @test length(out_basisset.basis[Li]) == 1
    end

    @testset "build_data dimensions match basisset" begin
        metadata = make_canonical_block_real_metadata(; hermitian=false)
        data = Quoll.build_data(metadata)

        basisset = Quoll.op_basisset(metadata)
        species2nbasis = Quoll.get_species2nbasis(basisset)

        body = Quoll.unwrap_data(data)
        for ((z1, z2), arr) in pairs(body)
            @test size(arr, 1) == species2nbasis[z1]
            @test size(arr, 2) == species2nbasis[z2]
        end
    end

end
