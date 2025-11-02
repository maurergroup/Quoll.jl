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
        @test sort(block_sparsity.images) == sort(ref_images)
        @test sort(block_sparsity.ij2images[(1, 1)]) == sort(ref_ij2images[(1, 1)])
        @test sort(block_sparsity.ij2images[(1, 2)]) == sort(ref_ij2images[(1, 2)])
        @test sort(block_sparsity.ij2images[(1, 3)]) == sort(ref_ij2images[(1, 3)])
        @test sort(block_sparsity.ij2images[(2, 1)]) == sort(ref_ij2images[(2, 1)])
        @test sort(block_sparsity.ij2images[(2, 2)]) == sort(ref_ij2images[(2, 2)])
        @test (2, 3) ∉ keys(block_sparsity.ij2images)
        @test sort(block_sparsity.ij2images[(3, 1)]) == sort(ref_ij2images[(3, 1)])
        @test (3, 2) ∉ keys(block_sparsity.ij2images)
        @test sort(block_sparsity.ij2images[(3, 3)]) == sort(ref_ij2images[(3, 3)])

        block_sparsity = Quoll.RealBlockSparsity(atoms, radii, hermitian = true)
        @test sort(block_sparsity.ij2images[(1, 1)]) == sort(ref_ij2images[(1, 1)])
        @test sort(block_sparsity.ij2images[(1, 2)]) == sort(ref_ij2images[(1, 2)])
        @test sort(block_sparsity.ij2images[(1, 3)]) == sort(ref_ij2images[(1, 3)])
        @test (2, 1) ∉ keys(block_sparsity.ij2images)
        @test sort(block_sparsity.ij2images[(2, 2)]) == sort(ref_ij2images[(2, 2)])
        @test (2, 3) ∉ keys(block_sparsity.ij2images)
        @test (3, 1) ∉ keys(block_sparsity.ij2images)
        @test (3, 2) ∉ keys(block_sparsity.ij2images)
        @test sort(block_sparsity.ij2images[(3, 3)]) == sort(ref_ij2images[(3, 3)])
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
        @test sort(block_sparsity.images) == sort(ref_images)
        @test sort(block_sparsity.ij2images[(1, 1)]) == sort(ref_ij2images[(1, 1)])
        @test sort(block_sparsity.ij2images[(1, 2)]) == sort(ref_ij2images[(1, 2)])
        @test (2, 1) ∉ keys(block_sparsity.ij2images)
        @test sort(block_sparsity.ij2images[(2, 2)]) == sort(ref_ij2images[(2, 2)])
    end

end

@testset "get_iexternal2ilocal" begin
    ij2images = dictionary([
        (1, 1) => [SA[0, 0, 1], SA[0, 0, -1], SA[0, 0, 0], SA[-1, 1, 0], SA[1, -1, 0]],
        (1, 2) => [SA[0, 0, 1], SA[0, 0, 0]],
        (2, 1) => [SA[0, 0, 0], SA[0, 0, -1]],
        (2, 2) => [SA[0, 0, 1], SA[0, 0, -1], SA[0, 0, 0], SA[-1, 1, 0], SA[1, -1, 0]],
    ])
    images = [
        SA[ 0,  0,  0],
        SA[ 0,  0,  1],
        SA[ 0,  0, -1],
        SA[-1,  1,  0],
        SA[ 1, -1,  0],
    ]
    hermitian = true
    sparsity = Quoll.RealBlockSparsity(ij2images, images, hermitian)
    n = nothing

    ref_iexternal2ilocal = dictionary([
        (1, 1) => [3, 1, 2, 4, 5],
        (1, 2) => [2, 1, n, n, n],
        (2, 1) => [1, n, 2, n, n],
        (2, 2) => [3, 1, 2, 4, 5],
    ])
    iexternal2ilocal = Quoll.get_iexternal2ilocal(images, sparsity)
    @test sort(keys(iexternal2ilocal)) == sort(keys(ref_iexternal2ilocal))
    for ij in keys(ref_iexternal2ilocal)
        @test iexternal2ilocal[ij] == ref_iexternal2ilocal[ij]
    end

    ref_iexternal2imlocal = dictionary([
        (1, 1) => [3, 2, 1, 5, 4],
        (1, 2) => [1, 2, n, n, n],
        (2, 1) => [2, n, 1, n ,n],
        (2, 2) => [3, 2, 1, 5, 4],
    ])

    iexternal2imlocal = Quoll.get_iexternal2imlocal(images, sparsity)
    @test sort(keys(iexternal2ilocal)) == sort(keys(ref_iexternal2ilocal))
    for ij in keys(ref_iexternal2ilocal)
        @test iexternal2ilocal[ij] == ref_iexternal2ilocal[ij]
    end
end

@testset "get_ilocal2imlocal" begin
    nhij2images = dictionary([
        (1, 1) => [SA[0, 0, 1], SA[0, 0, -1], SA[0, 0, 0], SA[-1, 1, 0], SA[1, -1, 0]],
        (1, 2) => [SA[0, 0, 1], SA[0, 0, 0]],
        (2, 1) => [SA[0, 0, 0], SA[0, 0, -1]],
        (2, 2) => [SA[0, 0, 1], SA[0, 0, 0], SA[-1, 1, 0], SA[1, -1, 0], SA[0, 0, -1]],
    ])
    images = [
        SA[ 0,  0,  0],
        SA[ 0,  0,  1],
        SA[ 0,  0, -1],
        SA[-1,  1,  0],
        SA[ 1, -1,  0],
    ]
    nhsparsity = Quoll.RealBlockSparsity(nhij2images, images, false)
    n = nothing

    @testset "Onsite" begin
        ref_onsite_ilocal2imlocal = dictionary([
            1 => [2, 1, 3, 5, 4],
            2 => [5, 2, 4, 3, 1],
        ])

        onsite_ilocal2imlocal = Quoll.get_onsite_ilocal2imlocal(nhsparsity)
        @test sort(keys(onsite_ilocal2imlocal)) == sort(keys(ref_onsite_ilocal2imlocal))
        for iat in keys(ref_onsite_ilocal2imlocal)
            @test onsite_ilocal2imlocal[iat] == ref_onsite_ilocal2imlocal[iat]
        end
    end

    @testset "General non-hermitian" begin
        ref_ilocal2imlocal = dictionary([
            (1, 1) => [2, 1, 3, 5, 4],
            (1, 2) => [2, 1],
            (2, 1) => [2, 1],
            (2, 2) => [5, 2, 4, 3, 1],
        ])

        ilocal2imlocal = Quoll.get_ilocal2imlocal(nhsparsity)
        @test sort(keys(ilocal2imlocal)) == sort(keys(ref_ilocal2imlocal))
        for iat in keys(ref_ilocal2imlocal)
            @test ilocal2imlocal[iat] == ref_ilocal2imlocal[iat]
        end
    end

end

@testset "convert_to_hermitian" begin
    ij2images = dictionary([
        (1, 1) => [SA[0, 0, 0], SA[0, 0, 1], SA[0, 1, 0], SA[0, 0, -1], SA[0, -1, 0]],
        (1, 2) => [SA[0, 0, 0], SA[0, 0, 1], SA[0, 1, 0]],
        (2, 1) => [SA[0, 0, 0], SA[0, 0, -1], SA[0, -1, 0]],
        (2, 2) => [SA[0, 0, 0], SA[0, 0, 1], SA[0, 1, 0], SA[0, 0, -1], SA[0, -1, 0]],
    ])
    images = [
        SA[0,  0,  0],
        SA[0,  0,  1],
        SA[0,  1,  0],
        SA[0,  0, -1],
        SA[0, -1,  0],
    ]
    hermitian = false
    nhsparsity = Quoll.RealBlockSparsity(ij2images, images, hermitian)
    hsparsity = Quoll.convert_to_hermitian(nhsparsity)

    ref_ij2images = dictionary([
        (1, 1) => [SA[0, 0, 0], SA[0, 0, 1], SA[0, 1, 0], SA[0, 0, -1], SA[0, -1, 0]],
        (1, 2) => [SA[0, 0, 0], SA[0, 0, 1], SA[0, 1, 0]],
        (2, 2) => [SA[0, 0, 0], SA[0, 0, 1], SA[0, 1, 0], SA[0, 0, -1], SA[0, -1, 0]],
    ])
    @test sort(keys(hsparsity.ij2images)) == sort(keys(ref_ij2images))
    @test sort(hsparsity.images) == sort(images)
    for ij in keys(ref_ij2images)
        @test sort(hsparsity.ij2images[ij]) == sort(ref_ij2images[ij])
    end
end

@testset "convert_to_nonhermitian" begin

    @testset "Missing pair" begin
        ij2images = dictionary([
            (1, 1) => [SA[0, 0, 0], SA[0, 0, 1], SA[0, 1, 0], SA[0, 0, -1], SA[0, -1, 0]],
            (1, 2) => [SA[0, 0, 0], SA[0, 0, 1], SA[0, 1, 0]],
            (2, 2) => [SA[0, 0, 0], SA[0, 0, 1], SA[0, 1, 0], SA[0, 0, -1], SA[0, -1, 0]],
        ])
        images = [
            SA[0,  0,  0],
            SA[0,  0,  1],
            SA[0,  1,  0],
            SA[0,  0, -1],
            SA[0, -1,  0],
        ]
        hermitian = true
        hsparsity = Quoll.RealBlockSparsity(ij2images, images, hermitian)
        nhsparsity = Quoll.convert_to_nonhermitian(hsparsity)

        ref_ij2images = dictionary([
            (1, 1) => [SA[0, 0, 0], SA[0, 0, 1], SA[0, 1, 0], SA[0, 0, -1], SA[0, -1, 0]],
            (1, 2) => [SA[0, 0, 0], SA[0, 0, 1], SA[0, 1, 0]],
            (2, 1) => [SA[0, 0, 0], SA[0, 0, -1], SA[0, -1, 0]],
            (2, 2) => [SA[0, 0, 0], SA[0, 0, 1], SA[0, 1, 0], SA[0, 0, -1], SA[0, -1, 0]],
        ])
        @test sort(keys(nhsparsity.ij2images)) == sort(keys(ref_ij2images))
        @test sort(nhsparsity.images) == sort(images)
        for ij in keys(ref_ij2images)
            @test sort(nhsparsity.ij2images[ij]) == sort(ref_ij2images[ij])
        end
    end

    @testset "Mixed" begin
        ij2images = dictionary([
            (1, 1) => [SA[0, 0, 0], SA[0, 0, 1], SA[0, 1, 0], SA[0, 0, -1], SA[0, -1, 0]],
            (1, 2) => [SA[0, 0, 0], SA[0, 0, 1]],
            (2, 1) => [SA[0, 0, 0], SA[0, -1, 0]],
            (2, 2) => [SA[0, 0, 0], SA[0, 0, 1], SA[0, 1, 0], SA[0, 0, -1], SA[0, -1, 0]],
        ])
        images = [
            SA[0,  0,  0],
            SA[0,  0,  1],
            SA[0,  1,  0],
            SA[0,  0, -1],
            SA[0, -1,  0],
        ]
        hermitian = true
        hsparsity = Quoll.RealBlockSparsity(ij2images, images, hermitian)
        nhsparsity = Quoll.convert_to_nonhermitian(hsparsity)

        # On-site still contain redundant images
        ref_ij2images = dictionary([
            (1, 1) => [SA[0, 0, 0], SA[0, 0, 1], SA[0, 1, 0], SA[0, 0, -1], SA[0, -1, 0]],
            (1, 2) => [SA[0, 0, 0], SA[0, 0, 1], SA[0, 1, 0]],
            (2, 1) => [SA[0, 0, 0], SA[0, 0, -1], SA[0, -1, 0]],
            (2, 2) => [SA[0, 0, 0], SA[0, 0, 1], SA[0, 1, 0], SA[0, 0, -1], SA[0, -1, 0]],
        ])
        @test sort(keys(nhsparsity.ij2images)) == sort(keys(ref_ij2images))
        @test sort(nhsparsity.images) == sort(images)
        for ij in keys(ref_ij2images)
            @test sort(nhsparsity.ij2images[ij]) == sort(ref_ij2images[ij])
        end
    end

end

@testset "make_onsite_hermitian!" begin
    ij2images = dictionary([
        (1, 1) => [SA[0, 0, 0], SA[1, 0, 0]],
        (1, 2) => [SA[0, 0, 0], SA[1, 0, 0]],
        (2, 2) => [SA[0, 0, 0], SA[0, 0, 1]],
    ])
    n_atoms = 2

    ij2images_ref = dictionary([
        (1, 1) => [SA[0, 0, 0], SA[1, 0, 0], SA[-1, 0, 0]],
        (1, 2) => [SA[0, 0, 0], SA[1, 0, 0]],
        (2, 2) => [SA[0, 0, 0], SA[0, 0, 1], SA[0, 0, -1]],
    ])

    Quoll.make_onsite_hermitian!(ij2images, n_atoms)
    @test sort(ij2images[(1, 1)]) == sort(ij2images_ref[(1, 1)])
    @test sort(ij2images[(2, 2)]) == sort(ij2images_ref[(2, 2)])
end