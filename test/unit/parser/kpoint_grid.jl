using Quoll
using Test
using Configurations
using StaticArrays

using Main.TestUtils

@testset "KPointGridParams" begin

    @testset "From parameters" begin
        d = Dict("mesh" => [10, 10, 1])
        params = from_dict(Quoll.Parser.KPointGridParams, d)
        @test params.kpoints === nothing
        @test params.mesh == SA[10, 10, 1]

        d = Dict("density" => 15.0)
        params = from_dict(Quoll.Parser.KPointGridParams, d)
        @test params.mesh === nothing
        @test params.kpoints === nothing
        @test params.density == 15.0

        d = Dict("shift" => true)
        params = from_dict(Quoll.Parser.KPointGridParams, d)
        @test params.shift == SA[true, true, true]

        d = Dict("shift" => [true, false, false])
        params = from_dict(Quoll.Parser.KPointGridParams, d)
        @test params.shift == SA[true, false, false]
    end

    @testset "From file" begin
        setupteardown_tmp() do

            kpoint_file = "kpoints.dat"
            open(kpoint_file, "w") do io
                write(io, "1.0 1.0 1.0 1.0\n")
                write(io, "1.0 0.0 1.0 0.0\n")
            end

            d = Dict("kpoints" => kpoint_file)
            params = from_dict(Quoll.Parser.KPointGridParams, d)
            @test params.kpoints == [SA[1.0, 1.0, 1.0, 1.0], SA[1.0, 0.0, 1.0, 0.0]]
        end
    end

    @testset "Invalid inputs" begin
        d = Dict("mesh" => [10, 10, -1])
        @test_throws ArgumentError from_dict(Quoll.Parser.KPointGridParams, d)

        d = Dict("mesh" => [10, 10])
        @test_throws DimensionMismatch from_dict(Quoll.Parser.KPointGridParams, d)

        d = Dict("density" => -10.0)
        @test_throws ArgumentError from_dict(Quoll.Parser.KPointGridParams, d)
    end

end