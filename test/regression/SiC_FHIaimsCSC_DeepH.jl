using DelimitedFiles
using JSON
using HDF5

using Test
using Main.TestUtils

np = 2
project = joinpath(@__DIR__, "../../app")
app = joinpath(@__DIR__, "../../app/quoll.jl")

inputfile = """
[input]
format = "FHI-aims"
directory = "SiC_FHIaims/SiC_FHIaims"
operators = ["H", "S"]

[output]
format = "DeepH"
directory = "SiC_DeepH_converted"
"""

SiC_examples = joinpath(@__DIR__, "../../examples/SiC")
tarballs = [
    joinpath(SiC_examples, "SiC_FHIaims.tar.gz")
    joinpath(SiC_examples, "SiC_DeepH.tar.gz")
]

setupteardown_tmp(tarballs = tarballs) do
    inputfile_path = "inputfile.toml"
    write(inputfile_path, inputfile)
    run.(mpiexec_quollapp(app, project, inputfile_path; np = np))

    dir_ref = joinpath("SiC_DeepH", "SiC_DeepH")
    dir = "SiC_DeepH_converted"

    @testset "hamiltonians.h5" begin
        db_ref = h5open(joinpath(dir_ref, "hamiltonians.h5"))
        db = h5open(joinpath(dir, "hamiltonians.h5"))

        k_ref = sort(keys(db_ref))
        k = sort(keys(db))

        blocks_ref = read.(Ref(db_ref), k_ref)
        blocks = read.(Ref(db), k)

        close(db_ref)
        close(db)

        @test k == k_ref
        @test blocks ≈ blocks_ref
    end

    @testset "overlaps.h5" begin
        db_ref = h5open(joinpath(dir_ref, "overlaps.h5"))
        db = h5open(joinpath(dir, "overlaps.h5"))

        k_ref = sort(keys(db_ref))
        k = sort(keys(db))

        blocks_ref = read.(Ref(db_ref), k_ref)
        blocks = read.(Ref(db), k)

        close(db_ref)
        close(db)

        @test k == k_ref
        @test blocks ≈ blocks_ref
    end

    # The way R_list is computed in DeepH-pack is different
    # (that method is not consistent with the images inside the `hamiltonians.h5`)
    # The method used in Quoll contains all the images (consistent with the true sparsity)
    @testset "R_list.dat" begin
        R_list_ref = readdlm(joinpath(dir_ref, "R_list.dat"))
        R_list = readdlm(joinpath(dir, "R_list.dat"))
        @test R_list ⊇ R_list_ref
    end

    @testset "element.dat" begin
        element_ref = readdlm(joinpath(dir_ref, "element.dat"))
        element = readdlm(joinpath(dir, "element.dat"))
        @test element == element_ref
    end

    @testset "orbital_types.dat" begin
        orbital_types_ref = readdlm(joinpath(dir_ref, "orbital_types.dat"))
        orbital_types = readdlm(joinpath(dir, "orbital_types.dat"))
        @test orbital_types == orbital_types_ref
    end

    @testset "lat.dat" begin
        lat_ref = readdlm(joinpath(dir_ref, "lat.dat"))
        lat = readdlm(joinpath(dir, "lat.dat"))
        @test lat ≈ lat_ref
    end

    @testset "rlat.dat" begin
        rlat_ref = readdlm(joinpath(dir_ref, "rlat.dat"))
        rlat = readdlm(joinpath(dir, "rlat.dat"))
        @test rlat ≈ rlat_ref
    end

    @testset "info.json" begin
        info_ref = JSON.parsefile(joinpath(dir_ref, "info.json"))
        info = JSON.parsefile(joinpath(dir, "info.json"))
        @test info == info_ref
    end

    @testset "site_positions.dat" begin
        site_positions_ref = readdlm(joinpath(dir_ref, "site_positions.dat"))
        site_positions = readdlm(joinpath(dir, "site_positions.dat"))
        @test get_translation_invariant(site_positions) ≈ get_translation_invariant(site_positions_ref)
    end

end