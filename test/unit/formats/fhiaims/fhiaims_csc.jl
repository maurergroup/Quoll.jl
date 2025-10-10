using Quoll
using StaticArrays
using DelimitedFiles
using HDF5

using Test
using Main.TestUtils

@testset "find_operatorkinds" begin

    paths = ["rs_hamiltonian.h5", "rs_overlap.out"]
    root, abspaths = create_temptree(paths)
    found_operatorkinds = find_operatorkinds(root, FHIaimsCSCOperator)
    @test Set(found_operatorkinds) == Set([Overlap(:ref), Hamiltonian(:ref)])

    paths = ["rs_overlap.out"]
    root, abspaths = create_temptree(paths)
    found_operatorkinds = find_operatorkinds(root, FHIaimsCSCOperator)
    @test Set(found_operatorkinds) == Set([Overlap(:ref)])

    paths = ["foo"]
    root, abspaths = create_temptree(paths)
    found_operatorkinds = find_operatorkinds(root, FHIaimsCSCOperator)
    @test isempty(found_operatorkinds)

end

@testset "RealCSCSparsity" begin
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
        sparsity = Quoll.RealCSCSparsity(pwd(), FHIaimsCSCOperator)
        @test sparsity.rowval == [1, 1, 2, 1, 1, 1]
        @test sparsity.colcellptr[1, :, :] == [1 2; 4 0; 5 6]
        @test sparsity.colcellptr[2, :, :] == [1 3; 4 -1; 5 6]
        @test sparsity.images == [SA[0, 0, 0], SA[0, 0, 1], SA[0, 0, -1]]
    end
end

@testset "load_operator_data" begin
    rs_hamiltonian = Float64[-1, -2, -3, -4, -4, -5, 0]

    @testset "Plain" begin
        setupteardown_tmp() do
            open("rs_hamiltonian.out", "w") do io
                writedlm(io, rs_hamiltonian)
            end
            operator_data = Quoll.load_operator_data(pwd(), Hamiltonian(:ref), FHIaimsCSCOperator)
            @test operator_data ≈ rs_hamiltonian[begin:end-1]
        end
    end

    @testset "HDF5" begin
        setupteardown_tmp() do
            h5open("rs_hamiltonian.h5", "w") do io
                dset = create_dataset(io, "sparse_matrix", Float64, (7,))
                write(dset, rs_hamiltonian)
            end
            operator_data = Quoll.load_operator_data(pwd(), Hamiltonian(:ref), FHIaimsCSCOperator)
            @test operator_data ≈ rs_hamiltonian[begin:end-1]
        end

    end

end
