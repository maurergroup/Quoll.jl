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
        # |x₂---||--x₁-|-----|| 
        cellvecs = SA[
            SA[1.0, 0.0, 0.0],
            SA[0.0, 100.0, 0.0],
            SA[0.0, 0.0, 100.0],
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
                ChemicalSpecies(:H1) => SA[0.40 - 0.00, 0.0, 0.0]u"Å",
                ChemicalSpecies(:H2) => SA[0.55 - 1.00, 0.0, 0.0]u"Å",
            ],
            cellvecs,
        )

        test_approx_eq(Quoll.recentre(atoms), recentered_atoms)
    end

    @testset "Non-periodic system" begin
        atoms = isolated_system([
            ChemicalSpecies(:H1) => SA[0.40, 0.00, 0.00]u"Å",
            ChemicalSpecies(:H2) => SA[0.55, 0.00, 0.00]u"Å",
        ])
        @test_throws ArgumentError Quoll.recentre(atoms)
    end
end

@testset "get_recip_lattice" begin
    real_lattice = collect(I(3))
    recip_lattice = collect(2π * I(3))
    @test Quoll.get_recip_lattice(real_lattice) == recip_lattice
end

@testset "get_frac_positions" begin
    cellvecs = SA[
        SA[10.0, 0.0, 0.0],
        SA[0.0, 10.0, 0.0],
        SA[0.0, 0.0, 10.0],
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

@testset "get_species2atom" begin
    cellvecs = SA[
        SA[10.0, 0.0, 0.0],
        SA[0.0, 10.0, 0.0],
        SA[0.0, 0.0, 10.0],
    ]u"Å"
    atoms = periodic_system(
        [
            ChemicalSpecies(:He) => SA[0.00, 1.00, 0.00]u"Å",
            ChemicalSpecies(:H) => SA[1.00, 0.00, 0.00]u"Å",
            ChemicalSpecies(:H) => SA[2.00, 0.00, 0.00]u"Å",
            ChemicalSpecies(:He) => SA[0.00, 2.00, 0.00]u"Å",
            ChemicalSpecies(:Li) => SA[0.00, 0.00, 1.00]u"Å",
            ChemicalSpecies(:He) => SA[0.00, 3.00, 0.00]u"Å",
        ],
        cellvecs,
    )

    species2atom_ref = Dictionary(
        [ChemicalSpecies(:H), ChemicalSpecies(:He), ChemicalSpecies(:Li)],
        [[2, 3], [1, 4, 6], [5]],
    )

    @test sort(Quoll.get_species2atom(atoms)) == sort(species2atom_ref)
end

@testset "get_images" begin
    cellvecs = SA[
        SA[1.0, 0.0,   0.0],
        SA[0.0, 1.0,   0.0],
        SA[0.0, 0.0, 100.0],
    ]u"Å"
    atoms = periodic_system(
        [
            ChemicalSpecies(:H) => SA[1.05, 0.50, 1.00]u"Å",
            ChemicalSpecies(:H) => SA[1.50, 1.05, 1.00]u"Å",
        ],
        cellvecs,
    )

    images_ref = [SA[1, 0, 0], SA[1, 1, 0]]
    @test Quoll.get_images(atoms) == images_ref
end

@testset "get_translation_invariant" begin
    coords = [
        0.0 1.0 2.0
        0.0 0.0 0.0
        0.0 0.0 0.0
    ]
    # Centroid is (1, 0, 0); subtracting it centres the coordinates.
    ref = [
        -1.0 0.0 1.0
        0.0 0.0 0.0
        0.0 0.0 0.0
    ]
    @test Quoll.get_translation_invariant(coords) ≈ ref
    # A rigid shift leaves the translation-invariant coordinates unchanged.
    shifted = coords .+ [5.0, -3.0, 2.0]
    @test Quoll.get_translation_invariant(shifted) ≈ Quoll.get_translation_invariant(coords)
end

@testset "same_atoms" begin
    cellvecs = SA[
        SA[10.0, 0.0, 0.0],
        SA[0.0, 10.0, 0.0],
        SA[0.0, 0.0, 10.0],
    ]u"Å"
    atoms = periodic_system(
        [
            ChemicalSpecies(:H) => SA[0.00, 0.00, 0.00]u"Å",
            ChemicalSpecies(:O) => SA[1.00, 0.00, 0.00]u"Å",
            ChemicalSpecies(:H) => SA[2.00, 0.00, 0.00]u"Å",
        ],
        cellvecs,
    )

    @testset "Identical system" begin
        @test Quoll.same_atoms(atoms, atoms)
    end

    @testset "Rigidly translated system" begin
        shift = SA[3.00, 1.00, -2.00]u"Å"
        atoms_shifted = periodic_system(
            [
                ChemicalSpecies(:H) => SA[0.00, 0.00, 0.00]u"Å" + shift,
                ChemicalSpecies(:O) => SA[1.00, 0.00, 0.00]u"Å" + shift,
                ChemicalSpecies(:H) => SA[2.00, 0.00, 0.00]u"Å" + shift,
            ],
            cellvecs,
        )
        @test Quoll.same_atoms(atoms, atoms_shifted)
    end

    @testset "Different species" begin
        atoms_species = periodic_system(
            [
                ChemicalSpecies(:H) => SA[0.00, 0.00, 0.00]u"Å",
                ChemicalSpecies(:N) => SA[1.00, 0.00, 0.00]u"Å",
                ChemicalSpecies(:H) => SA[2.00, 0.00, 0.00]u"Å",
            ],
            cellvecs,
        )
        @test !Quoll.same_atoms(atoms, atoms_species)
    end

    @testset "Different cell" begin
        cellvecs_other = SA[
            SA[12.0, 0.0, 0.0],
            SA[0.0, 10.0, 0.0],
            SA[0.0, 0.0, 10.0],
        ]u"Å"
        atoms_cell = periodic_system(
            [
                ChemicalSpecies(:H) => SA[0.00, 0.00, 0.00]u"Å",
                ChemicalSpecies(:O) => SA[1.00, 0.00, 0.00]u"Å",
                ChemicalSpecies(:H) => SA[2.00, 0.00, 0.00]u"Å",
            ],
            cellvecs_other,
        )
        @test !Quoll.same_atoms(atoms, atoms_cell)
    end

    @testset "Non-rigidly perturbed positions" begin
        atoms_perturbed = periodic_system(
            [
                ChemicalSpecies(:H) => SA[0.00, 0.00, 0.00]u"Å",
                ChemicalSpecies(:O) => SA[1.00, 0.50, 0.00]u"Å",
                ChemicalSpecies(:H) => SA[2.00, 0.00, 0.00]u"Å",
            ],
            cellvecs,
        )
        @test !Quoll.same_atoms(atoms, atoms_perturbed)
    end
end