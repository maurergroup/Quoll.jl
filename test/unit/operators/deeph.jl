using Quoll
using StaticArrays
using Dictionaries
using AtomsBase
using Unitful
using HDF5

using Test
using AtomsBaseTesting
using Main.TestUtils

function write_get_deeph_ref(::Type{M}) where {M<:Quoll.DeepHBlockRealMetadata}
    #! format: off
    ks = [
        ( 0,  0,  0, 1, 1),
        ( 1,  0,  0, 1, 1),
        (-1,  0,  0, 1, 1),
        ( 0,  0,  1, 1, 2),
        ( 0,  0, -1, 2, 1),
        ( 0,  0,  0, 2, 2),
        ( 1,  0,  0, 2, 2),
        (-1,  0,  0, 2, 2),
    ]
    #! format: on
    blocks = [reshape(fill(10.0 * i, 4) .+ collect(1:4), 2, 2) for i in axes(ks, 1)]
    if M <: Quoll.DeepHBlockSpinRealMetadata
        blocks = [Complex{Float64}.(block) for block in blocks]
    end
    data = Quoll.wrap_data(M, Dictionary{NTuple{5,Int},Matrix{Float64}}(ks, blocks))

    ij2images = dictionary([
        (1, 1) => [SA[0, 0, 0], SA[1, 0, 0], SA[-1, 0, 0]],
        (1, 2) => [SA[0, 0, 1]],
        (2, 1) => [SA[0, 0, -1]],
        (2, 2) => [SA[0, 0, 0], SA[1, 0, 0], SA[-1, 0, 0]],
    ])
    images = [
        SA[0, 0, 0],
        SA[0, 0, -1],
        SA[0, 0, 1],
        SA[-1, 0, 0],
        SA[1, 0, 0],
    ]
    sparsity = Quoll.BlockRealSparsity(ij2images, images, false)

    h5open("hamiltonians.h5", "w") do db
        for (k, block) in zip(ks, blocks)
            k_str = string(collect(k))
            write(db, k_str, permutedims(block, (2, 1)))
        end
    end

    return data, sparsity
end

@testset "load_atoms" begin
    element_dat = """
    1
    1
    """

    site_positions_dat = """
    0.0 1.0
    0.0 2.0
    0.0 3.0
    """

    lat_dat = """
    2.0 0.0 0.0
    0.0 4.0 0.0
    0.0 0.0 6.0
    """

    cell_vectors = [
        SA[2.0, 0.0, 0.0],
        SA[0.0, 4.0, 0.0],
        SA[0.0, 0.0, 6.0],
    ]u"Å"
    ref_atoms = periodic_system(
        [
            :H => [0.0, 0.0, 0.0]u"Å"
            :H => [1.0, 2.0, 3.0]u"Å"
        ],
        cell_vectors,
    )

    setupteardown_tmp() do
        open("element.dat", "w") do io
            write(io, element_dat)
        end
        open("site_positions.dat", "w") do io
            write(io, site_positions_dat)
        end
        open("lat.dat", "w") do io
            write(io, lat_dat)
        end
        dir = pwd()
        atoms = Quoll.load_atoms(Quoll.DeepHSource(), dir)
        test_approx_eq(atoms, ref_atoms)
    end
end

@testset "BasisSetMetadata" begin
    orbital_types_dat = """
    0 0 1
    0 1
    """
    setupteardown_tmp() do
        open("orbital_types.dat", "w") do io
            write(io, orbital_types_dat)
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

        basis_metadata = Quoll.BasisSetMetadata(Quoll.DeepHSource(), pwd(), atoms)

        # Note the reordering of p functions
        @test basis_metadata.basis[z1] ==
            [
            Quoll.BasisMetadata(z1, 1, 0, 0, nothing)
            Quoll.BasisMetadata(z1, 2, 0, 0, nothing)
            Quoll.BasisMetadata(z1, 1, 1, 1, nothing)
            Quoll.BasisMetadata(z1, 1, 1, -1, nothing)
            Quoll.BasisMetadata(z1, 1, 1, 0, nothing)
        ]
        @test basis_metadata.basis[z2] ==
            [
            Quoll.BasisMetadata(z2, 1, 0, 0, nothing)
            Quoll.BasisMetadata(z2, 1, 1, 1, nothing)
            Quoll.BasisMetadata(z2, 1, 1, -1, nothing)
            Quoll.BasisMetadata(z2, 1, 1, 0, nothing)
        ]
        @test basis_metadata.atom2species == [z1, z2]
    end
end

@testset "DeepHBlockRealMetadata" begin
    @testset "BlockRealSparsity" begin
        setupteardown_tmp() do
            _, ref_sparsity = write_get_deeph_ref(Quoll.DeepHBlockNoSpinRealMetadata)
            ref_images = ref_sparsity.images
            ref_ij2images = ref_sparsity.ij2images

            kind = Quoll.Hamiltonian(; source=:ref)
            sparsity = Quoll.BlockRealSparsity(
                Quoll.DeepHBlockNoSpinRealMetadata, pwd(), kind
            )
            images = sparsity.images
            ij2images = sparsity.ij2images

            @test sort(ref_images) == sort(images)
            for i in 1:2, j in 1:2
                @test sort(sparsity.ij2images[(i, j)]) == sort(ref_ij2images[(i, j)])
            end
        end
    end

    @testset "load_data" begin

        @testset "DeepHBlockNoSpinRealMetadata" begin
            setupteardown_tmp() do
                ref_data, _ = write_get_deeph_ref(Quoll.DeepHBlockNoSpinRealMetadata)
                ref_data_body = Quoll.unwrap_data(ref_data)

                kind = Quoll.Hamiltonian(; source=:ref)
                data = Quoll.load_data(Quoll.DeepHBlockNoSpinRealMetadata, pwd(), kind)
                data_body = Quoll.unwrap_data(data)
                
                @test keytype(data_body) <: NTuple{5,Int}
                @test eltype(data_body) <: Matrix{Float64}

                @test sort(keys(data_body)) == sort(keys(ref_data_body))
                for key in keys(ref_data_body)
                    @test data_body[key] == ref_data_body[key]
                end
            end
        end

        @testset "DeepHBlockSpinRealMetadata" begin
            setupteardown_tmp() do
                ref_data, _ = write_get_deeph_ref(Quoll.DeepHBlockSpinRealMetadata)
                ref_data_body = Quoll.unwrap_data(ref_data)

                kind = Quoll.Hamiltonian(; source=:ref, spin=:soc)
                data = Quoll.load_data(Quoll.DeepHBlockSpinRealMetadata, pwd(), kind)
                data_body = Quoll.unwrap_data(data)
                
                @test keytype(data_body) <: NTuple{5,Int}
                @test eltype(data_body) <: Matrix{Complex{Float64}}

                @test sort(keys(data_body)) == sort(keys(ref_data_body))
                for key in keys(ref_data_body)
                    @test data_body[key] == ref_data_body[key]
                end
            end
        end
    end
end
