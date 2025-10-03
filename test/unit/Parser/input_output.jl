using Quoll
using Test
using Configurations

include("../../testutils.jl")

@testset "InputParams" begin

    @testset "operators" begin
        directories, _ = create_tempdirs([])

        @testset "Nominal case" begin
            d = Dict("format" => "FHI-aims", "directories" => directories, "operators" => "H_ref")
            params = from_dict(Quoll.Parser.InputParams, d)

            @test params.operators == [Hamiltonian(:ref)]

            d = Dict("format" => "FHI-aims", "directories" => directories, "operators" => ["H_ref"])
            params = from_dict(Quoll.Parser.InputParams, d)

            @test params.operators == [Hamiltonian(:ref)]
        end

        @testset "Incorrect kind" begin
            d = Dict("format" => "FHI-aims", "directories" => directories, "operators" => ["foo", "bar"])

            @test_throws ArgumentError params = from_dict(Quoll.Parser.InputParams, d)
        end

    end

    @testset "format" begin
        directories, _ = create_tempdirs([])

        @testset "Nominal case" begin
            d = Dict("format" => "FHI-aims", "directories" => directories)
            params = from_dict(Quoll.Parser.InputParams, d)

            @test params.format == Quoll.OperatorIO.FHIaimsOperator
        end

        @testset "Incorrect format" begin
            d = Dict("format" => "foo", "directories" => directories)

            @test_throws ArgumentError from_dict(Quoll.Parser.InputParams, d)
        end
    end

    @testset "directories" begin

        @testset "Single directory" begin
            absroot, _ = create_tempdirs([])

            @testset "Absolute path" begin
                directories = absroot
                reference = [absroot]

                d = Dict("format" => "FHI-aims", "directories" => directories)
                params = from_dict(Quoll.Parser.InputParams, d)

                @test params.directories == reference
            end

            @testset "Relative path" begin
                directories = relpath(absroot)
                reference = [absroot]

                d = Dict("format" => "FHI-aims", "directories" => directories)
                params = from_dict(Quoll.Parser.InputParams, d)

                @test params.directories == reference
            end

        end

        @testset "Multiple directories" begin
            absroot, abspaths =
                create_tempdirs([joinpath("dir", "dir1"), joinpath("dir", "dir2")])

            @testset "Absolute paths" begin
                directories, reference = absroot, abspaths

                d = Dict("format" => "FHI-aims", "directories" => directories)
                params = from_dict(Quoll.Parser.InputParams, d)

                @test params.directories == reference
            end

            @testset "Relative paths" begin
                directories, reference = relpath(absroot), abspaths

                d = Dict("format" => "FHI-aims", "directories" => directories)
                params = from_dict(Quoll.Parser.InputParams, d)

                @test params.directories == reference
            end

        end

    end

end

@testset "OutputParams" begin

    @testset "format" begin
        directory = abspath("bar")

        @testset "Correct format" begin
            d = Dict("format" => "DeepH", "directory" => directory)
            params = from_dict(Quoll.Parser.OutputParams, d)

            @test params.format == Quoll.OperatorIO.DeepHOperator
        end

        @testset "Incorrect format" begin
            d = Dict("format" => "foo", "directory" => directory)

            @test_throws ArgumentError from_dict(Quoll.Parser.OutputParams, d)
        end

    end

    @testset "directory" begin

        @testset "Absolute directory" begin
            directory, reference = abspath("bar"), abspath("bar")

            d = Dict("format" => "DeepH", "directory" => directory)
            params = from_dict(Quoll.Parser.OutputParams, d)

            @test params.directory == reference
        end

        @testset "Relative directory" begin
            directory, reference = "bar", abspath("bar")

            d = Dict("format" => "DeepH", "directory" => directory)
            params = from_dict(Quoll.Parser.OutputParams, d)

            @test params.directory == reference
        end

    end

end