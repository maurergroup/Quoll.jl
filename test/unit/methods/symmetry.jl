using Quoll
using LinearAlgebra
using StaticArrays
using Unitful
using AtomsBase

using Test

@testset "get_symmetry_rotations" begin
    cellvecs = SA[
        SA[10.0,  0.0,  0.0],
        SA[ 0.0, 10.0,  0.0],
        SA[ 0.0,  0.0, 10.0],
    ]u"Å"
    atoms = periodic_system(
        [
            ChemicalSpecies(:H1) => SA[0.00, 0.00, 0.00]u"Å",
            ChemicalSpecies(:H2) => SA[1.00, 1.00, 1.00]u"Å",
        ],
        cellvecs,
    )
    
    spglib_cell = Quoll.get_spglib_cell(atoms)

    @testset "Nominal case" begin
        rotations = Quoll.get_symmetry_rotations(spglib_cell; crystal_symmetry = true)
        @test length(rotations) > 1
    end

    @testset "No crystal symmetry" begin
        rotations = Quoll.get_symmetry_rotations(spglib_cell; crystal_symmetry = false)
        @test length(rotations) == 1
        @test rotations[1] == SMatrix{3, 3}(I(3))
    end

end
