using Quoll
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
