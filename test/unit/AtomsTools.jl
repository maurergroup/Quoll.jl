using Quoll
using Test
using AtomsBase
using AtomsBaseTesting
using StaticArrays
using Unitful

include("../testutils.jl")

@testset "recentre" begin

    @testset "Periodic system" begin
        # ||: [ 0.0, 1.0)
        #  |: [-0.5, 0.5)
        # |-----||--x₁-|x₂---|| →
        # |x₂---||--x₁-|-----|| →
        # |-----||x₂---|--x₁-||
        cell = SA[
            SA[1.0,   0.0,   0.0],
            SA[0.0, 100.0,   0.0],
            SA[0.0,   0.0, 100.0],
        ]u"Å"
        atoms = periodic_system(
            [
                ChemicalSpecies(:H1) => SA[0.40, 0.00, 0.00]u"Å",
                ChemicalSpecies(:H2) => SA[0.55, 0.00, 0.00]u"Å",
            ],
            cell,
        )
        recentered_atoms = periodic_system(
            [
                ChemicalSpecies(:H1) => SA[0.40 + 0.50, 50.0, 50.0]u"Å",
                ChemicalSpecies(:H2) => SA[0.55 - 0.50, 50.0, 50.0]u"Å",
            ],
            cell,
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

@testset "load_atoms" begin

    @testset "FHIaimsOperator" begin
        # ||: [ 0.0, 1.0)
        #  |: [-0.5, 0.5)
        # |-----||--x₁-|x₂---|| →
        # |x₂---||--x₁-|-----|| →
        # |-----||x₂---|--x₁-||
    
        geometry = """
        lattice_vector 1.0   0.0   0.0
        lattice_vector 0.0 100.0   0.0
        lattice_vector 0.0   0.0 100.0
        atom 0.40 0.00 0.00 H
        atom 0.55 0.00 0.00 H
        """
        cell_vectors = SA[
            SA[1.0,   0.0,   0.0],
            SA[0.0, 100.0,   0.0],
            SA[0.0,   0.0, 100.0],
        ]u"Å"
        recentered_atoms = periodic_system(
            [
                ChemicalSpecies(:H) => SA[0.40 + 0.50, 50.0, 50.0]u"Å",
                ChemicalSpecies(:H) => SA[0.55 - 0.50, 50.0, 50.0]u"Å",
            ],
            cell_vectors,
        )

        tarballs = ()
        setupteardown_tmp(tarballs) do
            open("geometry.in", "w") do io
                write(io, geometry)
            end
            dir = pwd()
            atoms = load_atoms(dir, FHIaimsOperator)
            test_approx_eq(atoms, recentered_atoms)
        end

    end    

end