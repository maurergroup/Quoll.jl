using Quoll
using Test
using StaticArrays
using AtomsBase
using Unitful

include("../testutils.jl")

@testset "BasisSetMetadata" begin

    @testset "From file" begin
    
        @testset "FHIaims" begin
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

                n_basis = 2
                row_ind = [1, 1, 2, 1, 1, 1]
                col_cell_ptr = zeros(Int, 2, 3, 2)
                col_cell_ptr[1, :, :] = [1 2; 4 0; 5 6]
                col_cell_ptr[2, :, :] = [1 3; 4 -1; 5 6]
                cells = [SA[0, 0, 0], SA[0, 0, 1], SA[0, 0, -1]]

                operator_metadata = FHIaimsOperatorMetadata(n_basis, row_ind, col_cell_ptr, cells)
                basis_metadata = BasisSetMetadata(pwd(), atoms, operator_metadata)
                
                @test basis_metadata.basis[z1] == [BasisMetadata(z1, 1, 0, 0, Dict("type" => "atomic"))]
                @test basis_metadata.basis[z2] == [BasisMetadata(z2, 2, 1, -1, Dict("type" => "ionic"))]
                @test basis_metadata.n_basis_atom == [1, 1]
                @test basis_metadata.atom2species == [z1, z2]
                @test basis_metadata.basis2atom == [1, 2]
                @test basis_metadata.atom2basis == [1:1, 2:2]
            end
        end

    end

end