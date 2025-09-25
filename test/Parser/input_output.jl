using Quoll
using Test
using Configurations

include("../testutils.jl")

function test_from_dict(paramtype, directory, reference)
    d = Dict("format" => "foo", "directory" => directory)
    params = from_dict(paramtype, d)
    @test params.directory == reference
end

@testset "InputParams" begin

    @testset "Single directory" begin
        absroot, abspaths = create_tempdirs(["dir1"])

        @testset "Absolute path" begin
            test_from_dict(Quoll.Parser.InputParams, absroot, abspaths)
        end

        @testset "Relative path" begin
            test_from_dict(Quoll.Parser.InputParams, relpath(absroot), abspaths)
        end

    end

    @testset "Multiple directories" begin
        absroot, abspaths =
            create_tempdirs([joinpath("dir", "dir1"), joinpath("dir", "dir2")])

        @testset "Absolute paths" begin
            test_from_dict(Quoll.Parser.InputParams, absroot, abspaths)
        end

        @testset "Relative paths" begin
            test_from_dict(Quoll.Parser.InputParams, relpath(absroot), abspaths)
        end

    end

end

@testset "OutputParams" begin

    @testset "Absolute directory" begin
        test_from_dict(Quoll.Parser.OutputParams, abspath("bar"), abspath("bar"))
    end

    @testset "Relative directory" begin
        test_from_dict(Quoll.Parser.OutputParams, "bar", abspath("bar"))
    end

end