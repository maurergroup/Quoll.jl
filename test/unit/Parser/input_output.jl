using Quoll
using Test
using Configurations
using AtomsBase
using Unitful
using UnitfulAtomic

include("../../testutils.jl")

@testset "parse_radius" begin
    z = ChemicalSpecies(:Si)
    r = 10.0

    s = "Si 10"
    @test Quoll.Parser.parse_radius(s) == (z, r * u"Å")

    s = "Si 10.0"
    @test Quoll.Parser.parse_radius(s) == (z, r * u"Å")

    s = "Si 10.0 bohr"
    @test Quoll.Parser.parse_radius(s) == (z, uconvert(u"Å", r * u"bohr"))

    s = "Si 10.0 foo"
    @test_throws ArgumentError Quoll.Parser.parse_radius(s)

end


@testset "InputParams" begin

    @testset "operators" begin
        directories, _ = create_temptree([])

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
        directories, _ = create_temptree([])

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
            absroot, _ = create_temptree([])

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
                create_temptree([joinpath("dir", "dir1"), joinpath("dir", "dir2")])

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

    @testset "radii" begin
        directories, _ = create_temptree([])
    
        @testset "Nominal case" begin
            radii = ["Si 12.0 ang", "C 10.0 ang"]
            d = Dict("format" => "FHI-aims", "directories" => directories, "radii" => radii)
            params = from_dict(Quoll.Parser.InputParams, d)
            @test params.radii == Dict(
                ChemicalSpecies(:Si) => 12.0u"Å",
                ChemicalSpecies(:C) => 10.0u"Å",
            )
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