using Quoll
using AtomsBase
using Configurations

using Test
using Main.TestUtils

@testset "parse_basismetadata" begin
    @test collect(Quoll.Parser.parse_basismetadata(
        "C(1, 0, 0)"
    )) == [
        # BasisMetadata(ChemicalSpecies(:C), 1, 0, 0, Base.ImmutableDict{String, String}())
        BasisMetadata(ChemicalSpecies(:C), 1, 0, 0, Dict{String, String}())
    ]
    @test collect(Quoll.Parser.parse_basismetadata(
        "C(1, 0, 0) type = atomic"
    )) == [
        BasisMetadata(ChemicalSpecies(:C), 1, 0, 0, Dict("type" => "atomic"))
    ]
    @test collect(Quoll.Parser.parse_basismetadata(
        "Si(2, 1, *) type = atomic, func = gaussian"
    )) == [
        BasisMetadata(ChemicalSpecies(:Si), 2, 1, -1, Dict("type" => "atomic", "func" => "gaussian")),
        BasisMetadata(ChemicalSpecies(:Si), 2, 1,  0, Dict("type" => "atomic", "func" => "gaussian")),
        BasisMetadata(ChemicalSpecies(:Si), 2, 1,  1, Dict("type" => "atomic", "func" => "gaussian")),
    ]
end

@testset "BasisProjectionParams" begin
    projected_basis_str = [
        "C(1, 0, 0)  type = atomic",
        "Si(1, 0, 0) type = atomic",
        "Si(2, 1, *) type = hydrogenic, function = gaussian"
    ]
    projected_basis = [
        BasisMetadata(ChemicalSpecies(:C),  1, 0,  0, Dict("type" => "atomic")),
        BasisMetadata(ChemicalSpecies(:Si), 1, 0,  0, Dict("type" => "atomic")),
        BasisMetadata(ChemicalSpecies(:Si), 2, 1, -1, Dict("type" => "hydrogenic", "function" => "gaussian")),
        BasisMetadata(ChemicalSpecies(:Si), 2, 1,  0, Dict("type" => "hydrogenic", "function" => "gaussian")),
        BasisMetadata(ChemicalSpecies(:Si), 2, 1,  1, Dict("type" => "hydrogenic", "function" => "gaussian")),
    ]

    @testset "Existing method" begin
        d = Dict(
            "projected_basis" => projected_basis_str,
            "method" => "fc99v"
        )
        params = from_dict(Quoll.Parser.BasisProjectionParams, d)
        @test Set(params.projected_basis) == Set(projected_basis)
        @test params.method == FC99V
    end

    @testset "Non-existing method" begin
        d = Dict(
            "projected_basis" => projected_basis_str,
            "method" => "foo"
        )
        @test_throws ArgumentError from_dict(Quoll.Parser.BasisProjectionParams, d)
    end

end