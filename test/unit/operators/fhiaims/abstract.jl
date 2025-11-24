using Quoll
using StaticArrays
using AtomsBase
using Unitful
using UnitfulAtomic

using Test
using AtomsBaseTesting
using Main.TestUtils

@testset "load_atoms" begin

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
            Atom(:H, SA[0.40 + 0.50, 50.0, 50.0]u"Å", charge = 0u"e_au", magnetic_moment = 0),
            Atom(:H, SA[0.55 - 0.50, 50.0, 50.0]u"Å", charge = 0u"e_au", magnetic_moment = 0),
        ],
        cell_vectors,
    )

    setupteardown_tmp() do
        open("geometry.in", "w") do io
            write(io, geometry)
        end
        dir = pwd()
        atoms = load_atoms(dir, FHIaimsCSCOperator)
        test_approx_eq(atoms, recentered_atoms)
    end

end

@testset "BasisSetMetadata" begin
    basis_indices = """

    fn.   type   at.   n   l   m
        1 atomic     1   1   0   0
        2 ionic      2   2   1  -1
    """
    setupteardown_tmp() do
        open("basis-indices.out", "w") do io    
            write(io, basis_indices)
        end
        z1 = ChemicalSpecies(:H)
        z2 = ChemicalSpecies(:Li)
        cell_vectors = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]u"Å"
        atoms = periodic_system(
            [
                z1 => [0, 0, 0]u"Å"
                z2 => [0, 0, 0.5]u"Å"
            ],
            cell_vectors
        )

        basis_metadata = Quoll.BasisSetMetadata(pwd(), atoms, FHIaimsCSCOperator)
        
        @test basis_metadata.basis[z1] == [BasisMetadata(z1, 1, 0, 0, Dict("type" => "atomic"))]
        @test basis_metadata.basis[z2] == [BasisMetadata(z2, 2, 1, -1, Dict("type" => "ionic"))]
        @test basis_metadata.atom2species == [z1, z2]
    end
end
