using Quoll
using Configurations

using Test

@testset "SymmetryParams" begin

    @testset "Nominal case" begin
        d = Dict("space_group" => 21)
        params = from_dict(Quoll.Parser.SymmetryParams, d)
        @test params.space_group == 21
    end
    
    @testset "Invalid arguments" begin
        d = Dict("space_group" => 300)
        @test_throws ArgumentError from_dict(Quoll.Parser.SymmetryParams, d)
    end

end