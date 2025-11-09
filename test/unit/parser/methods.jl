using Quoll
using Configurations

using Test
using Main.TestUtils

struct DummyPostprocessParams
    fermi_level
    dos
end

struct DummyErrorParams
    mae
    eigenvalue_error
    el_entropy_error
end

struct DummyQuollParams
    basis_projection
    postprocessing
    error_metrics
end

# TODO: Test all the parameters that were read in once we know the file
# is unlikely to change anymore
# TODO: This would be more suitable as a regression test
# @testset "QuollParams" begin
#     SiC_exampledir = joinpath(@__DIR__, "../../../examples/SiC")
#     SiC_inputfile = joinpath(SiC_exampledir, "input_file.toml")
#     tarballs = [joinpath(SiC_exampledir, "SiC_FHIaims.tar.gz")]

#     setupteardown_tmp(tarballs = tarballs) do
#         params = from_toml(Quoll.Parser.QuollParams, SiC_inputfile)
#         @test params.input.format == FHIaimsCSCOperator
#     end
# end

@testset "get_output_dirs" begin

    @testset "Single root" begin
        @testset "Nominal case" begin
            output_root = abspath("output")
            input_roots = abspath.(["root"])
            input_eachroot = [
                abspath.([
                    joinpath("root", "dir1", "dir11"),
                    joinpath("root", "dir1", "dir12"),
                    joinpath("root", "dir2", "dir21"),
                    joinpath("root", "dir2", "dir22"),
                ]),
            ]

            output_dirs = Quoll.Parser.get_output_dirs(output_root, input_roots, input_eachroot)
            rel_output_dirs = relpath.(output_dirs, output_root)

            @test rel_output_dirs[1] == joinpath("dir1", "dir11")
            @test rel_output_dirs[2] == joinpath("dir1", "dir12")
            @test rel_output_dirs[3] == joinpath("dir2", "dir21")
            @test rel_output_dirs[4] == joinpath("dir2", "dir22")
        end

        @testset "Single directory" begin
            output_root = abspath("output")
            input_roots = abspath.(["root"])
            input_eachroot = [abspath.(["root"])]

            output_dirs = Quoll.Parser.get_output_dirs(output_root, input_roots, input_eachroot)
            rel_output_dirs = relpath.(output_dirs, output_root)

            @test rel_output_dirs[1] == "."
        end
    end

    @testset "Multiple roots" begin
        @testset "Nominal case" begin
            output_root = abspath("output")
            input_roots = abspath.(["root1", "root2"])
            input_eachroot = [
                abspath.([
                    joinpath("root1", "dir1", "dir11"),
                    joinpath("root1", "dir1", "dir12"),
                    joinpath("root1", "dir2", "dir21"),
                    joinpath("root1", "dir2", "dir22"),
                ]),
                abspath.([
                    joinpath("root2", "dir1", "dir11"),
                    joinpath("root2", "dir1", "dir12"),
                    joinpath("root2", "dir2", "dir21"),
                    joinpath("root2", "dir2", "dir22"),
                ]),
            ]

            output_dirs = Quoll.Parser.get_output_dirs(output_root, input_roots, input_eachroot)
            rel_output_dirs = relpath.(output_dirs, output_root)

            @test rel_output_dirs[1] == joinpath("root1", "dir1", "dir11")
            @test rel_output_dirs[2] == joinpath("root1", "dir1", "dir12")
            @test rel_output_dirs[3] == joinpath("root1", "dir2", "dir21")
            @test rel_output_dirs[4] == joinpath("root1", "dir2", "dir22")

            @test rel_output_dirs[5] == joinpath("root2", "dir1", "dir11")
            @test rel_output_dirs[6] == joinpath("root2", "dir1", "dir12")
            @test rel_output_dirs[7] == joinpath("root2", "dir2", "dir21")
            @test rel_output_dirs[8] == joinpath("root2", "dir2", "dir22")
        end

        @testset "Single directories" begin
            output_root = abspath("output")
            input_roots = abspath.(["root1", "root2"])
            input_eachroot = [abspath.(["root1"]), abspath.(["root2"])]

            output_dirs = Quoll.Parser.get_output_dirs(output_root, input_roots, input_eachroot)
            rel_output_dirs = relpath.(output_dirs, output_root)

            @test rel_output_dirs[1] == "root1"
            @test rel_output_dirs[2] == "root2"
        end
    end

end

@testset "find_leafdirs" begin

    @testset "Nominal case" begin
        tempdirs = [
            joinpath("dir1", "dir11"),
            joinpath("dir1", "dir12"),
            joinpath("dir2", "dir21"),
        ]

        root, paths = create_temptree(tempdirs)
        # root
        # ├── dir1
        # │   ├── dir11
        # │   └── dir12
        # └── dir2
        #     └── dir21

        leafdirs = Quoll.Parser.find_leafdirs(root)
        @test Set(leafdirs) == Set(paths)
    end

    @testset "Only root" begin
        tempdirs = []

        root, _ = create_temptree(tempdirs)
        # root
        # ├──

        leafdirs = Quoll.Parser.find_leafdirs(root)
        @test leafdirs == [root]

    end

end

@testset "requires_kpoint_grid" begin

    params = DummyQuollParams(
        nothing,
        DummyPostprocessParams(false, false),
        DummyErrorParams("foo", false, false),
    )
    @test Quoll.Parser.requires_kpoint_grid(params) == false
    
    params = DummyQuollParams(
        true,
        DummyPostprocessParams(false, false),
        DummyErrorParams("foo", false, false),
    )
    @test Quoll.Parser.requires_kpoint_grid(params) == true

    params = DummyQuollParams(
        nothing,
        DummyPostprocessParams(true, false),
        DummyErrorParams("foo", false, false),
    )
    @test Quoll.Parser.requires_kpoint_grid(params) == true

    params = DummyQuollParams(
        nothing,
        DummyPostprocessParams(false, false),
        DummyErrorParams("foo", true, false),
    )
    @test Quoll.Parser.requires_kpoint_grid(params) == true

    params = DummyQuollParams(
        true,
        DummyPostprocessParams(true, false),
        DummyErrorParams("foo", true, false),
    )
    @test Quoll.Parser.requires_kpoint_grid(params) == true
    
end

@testset "search_clashes" begin

    error_metrics = DummyErrorParams(false, false, false)
    @test isnothing(Quoll.Parser.search_clashes(nothing, error_metrics))

    error_metrics = DummyErrorParams(true, false, false)
    @test isnothing(Quoll.Parser.search_clashes(nothing, error_metrics))

    error_metrics = DummyErrorParams(false, false, false)
    @test isnothing(Quoll.Parser.search_clashes(true, error_metrics))

    error_metrics = DummyErrorParams(true, false, false)
    @test_throws ArgumentError Quoll.Parser.search_clashes(true, error_metrics)

end
