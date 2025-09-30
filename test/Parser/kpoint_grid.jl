using Quoll
using Test
using Configurations
using StaticArrays

@testset "KPointGridParams" begin

    @testset "Nominal case" begin
        d = Dict("grid" => [10, 10, 1])
        params = from_dict(Quoll.Parser.KPointGridParams, d)
        @test params.grid == SA[10, 10, 1]

        d = Dict("density" => 15.0)
        params = from_dict(Quoll.Parser.KPointGridParams, d)
        @test params.grid === nothing
        @test params.density == 15.0
    end

    @testset "Invalid inputs" begin
        d = Dict("grid" => [10, 10, -1])
        @test_throws ArgumentError from_dict(Quoll.Parser.KPointGridParams, d)

        d = Dict("grid" => [10, 10])
        @test_throws DimensionMismatch from_dict(Quoll.Parser.KPointGridParams, d)

        d = Dict("density" => -10.0)
        @test_throws ArgumentError from_dict(Quoll.Parser.KPointGridParams, d)
    end

end