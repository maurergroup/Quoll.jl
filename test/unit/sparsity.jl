using Quoll
using StaticArrays
using Dictionaries
using AtomsBase
using Unitful

using Test

@testset "RealBlockSparsity" begin

    @testset "From neighbourlist" begin
        # |H₁-H₂-----H₃-|H₁-H₂-----H₃-|H₁-H₂-----H₃-|
        # cutoff: --
        cellvecs = SA[
            SA[1.0,   0.0,   0.0],
            SA[0.0, 100.0,   0.0],
            SA[0.0,   0.0, 100.0],
        ]u"Å"
        atoms = periodic_system(
            [
                ChemicalSpecies(:H) => SA[0.05, 1.00, 1.00]u"Å",
                ChemicalSpecies(:H) => SA[0.15, 1.00, 1.00]u"Å",
                ChemicalSpecies(:H) => SA[0.85, 1.00, 1.00]u"Å",
            ],
            cellvecs,
        )
        radii = Dict(ChemicalSpecies(:H) => 0.11u"Å")

        ref_images = [SA[0, 0, 0], SA[-1, 0, 0], SA[1, 0, 0]]
        ref_ij2images = dictionary([
            (1, 1) => [SA[0, 0, 0]                           ],
            (1, 2) => [SA[0, 0, 0]                           ],
            (1, 3) => [             SA[-1, 0, 0]             ],
            (2, 1) => [SA[0, 0, 0]                           ],
            (2, 2) => [SA[0, 0, 0]                           ],
           #(2, 3) => [],
            (3, 1) => [                           SA[1, 0, 0]],
           #(3, 2) => [],
            (3, 3) => [SA[0, 0, 0]                           ],
        ])

        block_sparsity = Quoll.RealBlockSparsity(atoms, radii, hermitian = false)
        @test Set(block_sparsity.images) == Set(ref_images)
        @test Set(block_sparsity.ij2images[(1, 1)]) == Set(ref_ij2images[(1, 1)])
        @test Set(block_sparsity.ij2images[(1, 2)]) == Set(ref_ij2images[(1, 2)])
        @test Set(block_sparsity.ij2images[(1, 3)]) == Set(ref_ij2images[(1, 3)])
        @test Set(block_sparsity.ij2images[(2, 1)]) == Set(ref_ij2images[(2, 1)])
        @test Set(block_sparsity.ij2images[(2, 2)]) == Set(ref_ij2images[(2, 2)])
        @test (2, 3) ∉ keys(block_sparsity.ij2images)
        @test Set(block_sparsity.ij2images[(3, 1)]) == Set(ref_ij2images[(3, 1)])
        @test (3, 2) ∉ keys(block_sparsity.ij2images)
        @test Set(block_sparsity.ij2images[(3, 3)]) == Set(ref_ij2images[(3, 3)])

        block_sparsity = Quoll.RealBlockSparsity(atoms, radii, hermitian = true)
        @test Set(block_sparsity.ij2images[(1, 1)]) == Set(ref_ij2images[(1, 1)])
        @test Set(block_sparsity.ij2images[(1, 2)]) == Set(ref_ij2images[(1, 2)])
        @test Set(block_sparsity.ij2images[(1, 3)]) == Set(ref_ij2images[(1, 3)])
        @test (2, 1) ∉ keys(block_sparsity.ij2images)
        @test Set(block_sparsity.ij2images[(2, 2)]) == Set(ref_ij2images[(2, 2)])
        @test (2, 3) ∉ keys(block_sparsity.ij2images)
        @test (3, 1) ∉ keys(block_sparsity.ij2images)
        @test (3, 2) ∉ keys(block_sparsity.ij2images)
        @test Set(block_sparsity.ij2images[(3, 3)]) == Set(ref_ij2images[(3, 3)])
    end

    @testset "From RealCSCSparsity" begin
        rowval = [1, 1, 2, 1, 1, 1]
        colcellptr = zeros(Int, 2, 3, 2)
        colcellptr[1, :, :] = [1 2; 4 0; 5 6]
        colcellptr[2, :, :] = [1 3; 4 -1; 5 6]
        images = [SA[0, 0, 0], SA[0, 0, 1], SA[0, 0, -1]]
        basis2atom = [1, 2]

        ref_images = [SA[0, 0, 0], SA[0, 0, 1], SA[0, 0, -1]]
        ref_ij2images = dictionary([
            (1, 1) => [SA[0, 0, 0], SA[0, 0, 1], SA[0, 0, -1]]
            (1, 2) => [SA[0, 0, 0],              SA[0, 0, -1]]
           #(2, 1) => [SA[0, 0, 0], SA[0, 0, 1]              ]
            (2, 2) => [SA[0, 0, 0]                           ]
        ])

        block_sparsity = Quoll.RealBlockSparsity(colcellptr, rowval, images, basis2atom)
        @test block_sparsity.hermitian == true
        @test Set(block_sparsity.images) == Set(ref_images)
        @test Set(block_sparsity.ij2images[(1, 1)]) == Set(ref_ij2images[(1, 1)])
        @test Set(block_sparsity.ij2images[(1, 2)]) == Set(ref_ij2images[(1, 2)])
        @test (2, 1) ∉ keys(block_sparsity.ij2images)
        @test Set(block_sparsity.ij2images[(2, 2)]) == Set(ref_ij2images[(2, 2)])
    end

end