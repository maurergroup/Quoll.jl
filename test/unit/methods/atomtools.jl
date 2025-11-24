using Quoll
using LinearAlgebra
using Dictionaries
using StaticArrays
using AtomsBase
using Unitful

using Test
using AtomsBaseTesting
using Main.TestUtils

@testset "recentre" begin

    @testset "Periodic system" begin
        # ||: [ 0.0, 1.0)
        #  |: [-0.5, 0.5)
        # |-----||--x₁-|x₂---|| →
        # |x₂---||--x₁-|-----|| →
        # |-----||x₂---|--x₁-||
        cellvecs = SA[
            SA[1.0,   0.0,   0.0],
            SA[0.0, 100.0,   0.0],
            SA[0.0,   0.0, 100.0],
        ]u"Å"
        atoms = periodic_system(
            [
                ChemicalSpecies(:H1) => SA[0.40, 0.00, 0.00]u"Å",
                ChemicalSpecies(:H2) => SA[0.55, 0.00, 0.00]u"Å",
            ],
            cellvecs,
        )
        recentered_atoms = periodic_system(
            [
                ChemicalSpecies(:H1) => SA[0.40 + 0.50, 50.0, 50.0]u"Å",
                ChemicalSpecies(:H2) => SA[0.55 - 0.50, 50.0, 50.0]u"Å",
            ],
            cellvecs,
        )

        test_approx_eq(recentre(atoms), recentered_atoms)
    end

    @testset "Non-periodic system" begin
        atoms = isolated_system([
            ChemicalSpecies(:H1) => SA[0.40, 0.00, 0.00]u"Å",
            ChemicalSpecies(:H2) => SA[0.55, 0.00, 0.00]u"Å",
        ])
        @test_throws ArgumentError recentre(atoms)
    end

end

@testset "get_recip_lattice" begin
    real_lattice = collect(I(3))
    recip_lattice = collect(2π * I(3))
    @test Quoll.get_recip_lattice(real_lattice) == recip_lattice
end

@testset "get_frac_positions" begin
    cellvecs = SA[
        SA[10.0,  0.0,  0.0],
        SA[ 0.0, 10.0,  0.0],
        SA[ 0.0,  0.0, 10.0],
    ]u"Å"
    atoms = periodic_system(
        [
            ChemicalSpecies(:H1) => SA[0.50, 0.00, 0.00]u"Å",
            ChemicalSpecies(:H2) => SA[2.00, 0.00, 0.00]u"Å",
        ],
        cellvecs,
    )

    frac_positions = [0.05 0.20; 0.00 0.00; 0.00 0.00]
    @test Quoll.get_frac_positions(atoms) ≈ frac_positions
end
