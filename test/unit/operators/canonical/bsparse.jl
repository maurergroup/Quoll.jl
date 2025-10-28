using Quoll
using AtomsBase

using Test

# function get_z1z2_ij2interval(atoms::AbstractSystem, sparsity::RealBlockSparsity)
@testset "get_z1z2_ij2interval" begin
    cellvecs = SA[
        SA[1.0,   0.0,   0.0],
        SA[0.0, 100.0,   0.0],
        SA[0.0,   0.0, 100.0],
    ]u"Å"
    a, b = ChemicalSpecies(:H1), ChemicalSpecies(:H2)
    atoms = periodic_system(
        [
            a => SA[0.40, 0.00, 0.00]u"Å",
            a => SA[0.45, 0.00, 0.00]u"Å",
            b => SA[0.55, 0.00, 0.00]u"Å",
            b => SA[0.60, 0.00, 0.00]u"Å",
        ],
        cellvecs,
    )

    ij2images = dictionary([
        (1, 1) => [SA[0, 0, 0], SA[1, 0, 0], SA[-1, 0, 0]] # a, a
        (1, 2) => [SA[0, 0, 0], SA[1, 0, 0]]               # a, a
        (1, 3) => [SA[0, 0, 0], SA[1, 0, 0]]               # a, b
        (2, 2) => [SA[0, 0, 0], SA[0, 0, 1], SA[0, 0, -1]] # a, a
        (2, 3) => [SA[0, 0, 0], SA[0, 0, 1], SA[0, 0, -1]] # a, b
        (2, 4) => [SA[0, 0, 0], SA[0, 0, -1]]              # a, b
        (3, 4) => [SA[0, 0, 0], SA[0, 0, 1], SA[0, 0, -1]] # b, b
        (4, 4) => [SA[0, 0, 0], SA[0, 0, -1]]              # b, b
    ])
    images = [
        SA[ 0,  0,  0],
        SA[ 1,  0,  0],
        SA[ 0,  0,  1],
        SA[-1,  0,  0],
        SA[ 0,  0, -1],
    ]
    hermitian = true
    sparsity = Quoll.RealBlockSparsity(ij2images, images, hermitian)

    ref_z1z2_ij2interval = dictionary([
        (a, a) => dictionary([
            (1, 1) => 1:3
            (1, 2) => 4:5
            (2, 2) => 6:8
        ]),
        (a, b) => dictionary([
            (1, 3) => 1:2
            (2, 3) => 3:5
            (2, 4) => 6:7
        ]),
        (b, b) => dictionary([
            (3, 4) => 1:3
            (4, 4) => 4:5
        ])
    ])

    z1z2_ij2interval = Quoll.get_z1z2_ij2interval(atoms, sparsity)
    @test keys(z1z2_ij2interval) == keys(ref_z1z2_ij2interval)
    for z1z2 in keys(ref_z1z2_ij2interval)
        @test keys(z1z2_ij2interval[z1z2]) == keys(ref_z1z2_ij2interval[z1z2])
        for ij in keys(ref_z1z2_ij2interval[z1z2])
            @test z1z2_ij2interval[z1z2][ij] == ref_z1z2_ij2interval[z1z2][ij]
        end
    end
end
