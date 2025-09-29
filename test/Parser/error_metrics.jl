using Quoll
using Configurations
using Test

@testset "ElEntropyErrorParams" begin

    @testset "Nominal case" begin
        d = Dict(
            "smearing_function" => "gaussian",
            "temperature" => 500
        )
        params = from_dict(Quoll.Parser.ElEntropyErrorParams, d)
        @test params.smearing_function == Gaussian()
        @test params.temperature == 500.0
    end

    @testset "Unphysical temperature" begin
        d = Dict(
            "smearing_function" => "fermi_dirac",
            "temperature" => -100
        )
        @test_throws ArgumentError from_dict(Quoll.Parser.ElEntropyErrorParams, d)
    end

    @testset "Non-existing smearing" begin
        d = Dict(
            "smearing_function" => "foo",
            "temperature" => 1000
        )
        @test_throws ArgumentError from_dict(Quoll.Parser.ElEntropyErrorParams, d)
    end

end

@testset "EigenvalueErrorParams" begin

    d = Dict(
        "smearing_function" => "gaussian",
        "temperature" => 500
    )
    params = from_dict(Quoll.Parser.EigenvalueErrorParams, d)
    @test params.smearing_function == Gaussian()
    @test params.temperature == 500.0

end

@testset "ErrorMetricParams" begin

    @testset "Nominal case" begin
        d = Dict(
            "mae" => true,
            "eigenvalue_error" => true,
            "el_entropy_error" => true,
            "eigenvalue_error_params" => Dict(
                "smearing_function" => "gaussian",
                "temperature" => 500
            ),
            "el_entropy_error_params" => Dict(
                "smearing_function" => "gaussian",
                "temperature" => 500
            )
        )
        params = from_dict(Quoll.Parser.ErrorMetricParams, d)
        @test params.mae == true
        @test params.eigenvalue_error == true
        @test params.el_entropy_error == true
        @test params.eigenvalue_error_params.smearing_function == Gaussian()
        @test params.eigenvalue_error_params.temperature == 500.0
        @test params.el_entropy_error_params.smearing_function == Gaussian()
        @test params.el_entropy_error_params.temperature == 500.0
    end

    @testset "Default values" begin
        d = Dict{String, Any}()
        default_params = Quoll.Parser.ErrorMetricParams()
        params = from_dict(Quoll.Parser.ErrorMetricParams, d)
        @test params.mae                                       == default_params.mae                                      
        @test params.eigenvalue_error                          == default_params.eigenvalue_error                         
        @test params.el_entropy_error                          == default_params.el_entropy_error                         
        @test params.eigenvalue_error_params.smearing_function == default_params.eigenvalue_error_params.smearing_function
        @test params.eigenvalue_error_params.temperature       == default_params.eigenvalue_error_params.temperature      
        @test params.el_entropy_error_params.smearing_function == default_params.el_entropy_error_params.smearing_function
        @test params.el_entropy_error_params.temperature       == default_params.el_entropy_error_params.temperature      
    end

end