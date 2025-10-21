using Quoll

using Test

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