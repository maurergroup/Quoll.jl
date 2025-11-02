using Quoll
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

@testset "speciesdict_to_atomarray" begin
    z1 = ChemicalSpecies(:H1)
    z2 = ChemicalSpecies(:H2)
    atom2species = [z2, z1, z2, z1]

    vz1 = [5, 6]
    vz2 = [7, 8]

    d = dictionary([
        z1 => vz1,
        z2 => vz2,
    ])
    ref = [vz2, vz1, vz2, vz1]

    @test isequal(Quoll.speciesdict_to_atomarray(d, atom2species, DIM = 1), ref)

    vz1z1 = [10, 20]
    vz1z2 = [30, 40]
    vz2z1 = missing
    vz2z2 = [50, 60, 70]

    d = dictionary([
        (z1, z1) => vz1z1,
        (z1, z2) => vz1z2,
        (z2, z2) => vz2z2,
    ])
    ref = reshape([
        vz2z2, vz1z2, vz2z2, vz1z2,
        vz2z1, vz1z1, vz2z1, vz1z1,
        vz2z2, vz1z2, vz2z2, vz1z2,
        vz2z1, vz1z1, vz2z1, vz1z1,
    ], 4, 4)

    @test isequal(Quoll.speciesdict_to_atomarray(d, atom2species, DIM = 2), ref)
end