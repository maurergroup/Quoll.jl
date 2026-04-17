using Quoll
using StaticArrays
using AtomsBase
using Unitful
using UnitfulAtomic
using DelimitedFiles
using HDF5

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
        SA[1.0, 0.0, 0.0],
        SA[0.0, 100.0, 0.0],
        SA[0.0, 0.0, 100.0],
    ]u"Å"
    recentered_atoms = periodic_system(
        [
            Atom(:H, SA[0.40 + 0.50, 50.0, 50.0]u"Å"; charge=0u"e_au", magnetic_moment=0),
            Atom(:H, SA[0.55 - 0.50, 50.0, 50.0]u"Å"; charge=0u"e_au", magnetic_moment=0),
        ],
        cell_vectors,
    )

    setupteardown_tmp() do
        open("geometry.in", "w") do io
            write(io, geometry)
        end
        dir = pwd()
        atoms = Quoll.load_atoms(Quoll.FHIaimsSource(), dir)
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
            cell_vectors,
        )

        basis_metadata = Quoll.BasisSetMetadata(Quoll.FHIaimsSource(), pwd(), atoms)

        @test basis_metadata.basis[z1] ==
            [Quoll.BasisMetadata(z1, 1, 0, 0, Dict(:type => :atomic))]
        @test basis_metadata.basis[z2] ==
            [Quoll.BasisMetadata(z2, 2, 1, -1, Dict(:type => :ionic))]
        @test basis_metadata.atom2species == [z1, z2]
    end
end

@testset "FHIaimsCSCRealMetadata" begin
    @testset "find_operatorkinds" begin
        paths = ["rs_hamiltonian.h5", "rs_overlap.out"]
        root, abspaths = create_temptree(paths)
        found_operatorkinds = find_operatorkinds(Quoll.FHIaimsCSCNoSpinRealMetadata, root)
        @test Set(found_operatorkinds) ==
            Set([Overlap(; source=:ref), Hamiltonian(; source=:ref)])

        paths = ["rs_hamiltonian_up.h5", "rs_hamiltonian_down.out", "rs_overlap.h5"]
        root, abspaths = create_temptree(paths)
        found_operatorkinds = find_operatorkinds(Quoll.FHIaimsCSCSpinRealMetadata, root)
        @test Set(found_operatorkinds) ==
            Set([
            Overlap(; source=:ref, spin=:up),
            Overlap(; source=:ref, spin=:down),
            Hamiltonian(; source=:ref, spin=:up),
        ])

        paths = ["rs_overlap.out"]
        root, abspaths = create_temptree(paths)
        found_operatorkinds = find_operatorkinds(Quoll.FHIaimsCSCNoSpinRealMetadata, root)
        @test Set(found_operatorkinds) == Set([Overlap(; source=:ref)])

        paths = ["foo"]
        root, abspaths = create_temptree(paths)
        found_operatorkinds = find_operatorkinds(Quoll.FHIaimsCSCNoSpinRealMetadata, root)
        @test isempty(found_operatorkinds)
    end

    @testset "CSCRealSparsity" begin

        @testset "Periodic system" begin
            rs_indices = """
            n_hamiltonian_matrix_size: 7
            n_cells_in_hamiltonian: 4
            n_basis: 2
            cell_index
                    0          0          0
                    0          0          1
                    0          0         -1
            999999999  999999999  999999999
            index_hamiltonian(1,:,:)
                    1          2
                    4          0
                    5          6
                    0         -1
            index_hamiltonian(2,:,:)
                    1          3
                    4         -1
                    5          6
                    0         -1
            column_index_hamiltonian
                    1
                    1
                    2
                    1
                    1
                    1
                    0
            """
            setupteardown_tmp() do
                open("rs_indices.out", "w") do io
                    write(io, rs_indices)
                end
                sparsity = Quoll.CSCRealSparsity(Quoll.FHIaimsSource(), pwd())
                @test sparsity.rowval == [1, 1, 2, 1, 1, 1]
                @test sparsity.colcellptr[1, :, :] == [1 2; 4 0; 5 6]
                @test sparsity.colcellptr[2, :, :] == [1 3; 4 -1; 5 6]
                @test sparsity.images == [SA[0, 0, 0], SA[0, 0, 1], SA[0, 0, -1]]
            end
        end

        @testset "Isolated system" begin
            rs_indices = """
            n_hamiltonian_matrix_size: 4
            n_cells_in_hamiltonian: 2
            n_basis: 2
            cell_index (not used in aperiodic systems)
            index_hamiltonian(1,:,:)
                    1          2
                    0         -1
            index_hamiltonian(2,:,:)
                    1          3
                    0         -1
            column_index_hamiltonian
                    1
                    1
                    2
                    2
                    0
            """
            setupteardown_tmp() do
                open("rs_indices.out", "w") do io
                    write(io, rs_indices)
                end
                sparsity = Quoll.CSCRealSparsity(Quoll.FHIaimsSource(), pwd())
                @test sparsity.rowval == [1, 1, 2, 2]
                @test sparsity.colcellptr[1, :, :] == [1 2;]
                @test sparsity.colcellptr[2, :, :] == [1 3;]
                @test sparsity.images == [SA[0, 0, 0]]
            end
        end
    end

    @testset "load_data" begin
        rs_hamiltonian = Float64[-1, -2, -3, -4, -4, -5, 0]

        @testset "Plain" begin
            setupteardown_tmp() do
                open("rs_hamiltonian.out", "w") do io
                    writedlm(io, rs_hamiltonian)
                end
                operator_data = Quoll.load_data(
                    Quoll.FHIaimsCSCNoSpinRealMetadata, pwd(), Hamiltonian(; source=:ref)
                )
                @test Quoll.unwrap_data(operator_data) ≈
                    rs_hamiltonian[begin:(end - 1)] * Quoll.HARTREE_CODATA_2002
            end
        end

        @testset "HDF5" begin
            setupteardown_tmp() do
                h5open("rs_hamiltonian.h5", "w") do io
                    dset = create_dataset(io, "sparse_matrix", Float64, (7,))
                    write(dset, rs_hamiltonian)
                end
                operator_data = Quoll.load_data(
                    Quoll.FHIaimsCSCNoSpinRealMetadata, pwd(), Hamiltonian(; source=:ref)
                )
                @test Quoll.unwrap_data(operator_data) ≈
                    rs_hamiltonian[begin:(end - 1)] * Quoll.HARTREE_CODATA_2002
            end
        end
    end
end