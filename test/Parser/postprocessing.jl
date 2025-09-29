using Quoll
using Test
using Configurations
using StaticArrays

@testset "parse_kpathparams" begin

    @testset "Nominal case" begin
        band = "0.0  0.0  0.0 -0.5 -0.5 -0.5 20 G X"
        kpathparams = Quoll.Parser.parse_kpathparams(band)
        @test kpathparams.kpoint_begin == SA[0.0, 0.0, 0.0]
        @test kpathparams.kpoint_end == SA[-0.5, -0.5, -0.5]
        @test kpathparams.n_points == 20
        @test kpathparams.symbol_begin == :G
        @test kpathparams.symbol_end == :X
    end

    @testset "Invalid bands" begin
        band = "foo  0.0  0.0 -0.5 -0.5 -0.5 20 G X"
        @test_throws ArgumentError Quoll.Parser.parse_kpathparams(band)

        band = "0.0  0.0  0.0 -0.5 -0.5 -0.5 20.1 G X"
        @test_throws ArgumentError Quoll.Parser.parse_kpathparams(band)
        
        band = "0.0  0.0  0.0 -0.5 -0.5 -0.5 20 G"
        @test_throws ArgumentError Quoll.Parser.parse_kpathparams(band)
    end

end

@testset "DOSParams" begin

    @testset "Inverted energy ranges" begin
        d = Dict(
            "smearing_function" => "fermi_dirac",
            "temperature" => 500.0,
            "n_points" => 1000,
            "energy_begin" => 20.0,
            "energy_end" => -20.0,
        )
        @test_throws ArgumentError from_dict(Quoll.Parser.DOSParams, d)
    end

end

@testset "PostprocessParams" begin

    @testset "Nominal case" begin
        d = Dict(
            "fermi_level" => true,
            "dos" => false,
            "band_structure" => false,
            "fermi_level_params" => Dict(
                "smearing_function" => "fermi_dirac"
            ),
            "dos_params" => Dict(
                "smearing_function" => "fermi_dirac",
                "temperature" => 500.0,
                "n_points" => 1000,
                "energy_begin" => -20.0,
                "energy_end" => 20.0,
            ),
            "band_structure_params" => Dict(
                "bands" => [
                    "0.0  0.0  0.0 -0.5 -0.5 -0.5 20 G X",
                    "0.5  0.0  0.0  0.0  0.0  0.0 20 Y G"
                ]
            )
        )
        params = from_dict(Quoll.Parser.PostprocessParams, d)

        @test params.fermi_level                           == true
        @test params.dos                                   == false
        @test params.band_structure                        == false
        @test params.fermi_level_params.smearing_function  == FermiDirac()
        @test params.dos_params.smearing_function          == FermiDirac()
        @test params.dos_params.temperature                == 500.0
        @test params.dos_params.n_points                   == 1000
        @test params.dos_params.energy_begin               == -20.0
        @test params.dos_params.energy_end                 == 20.0
        @test params.band_structure_params.bands           == [
            Quoll.Parser.KPathParams(SA[0.0, 0.0, 0.0], SA[-0.5, -0.5, -0.5], 20, :G, :X),
            Quoll.Parser.KPathParams(SA[0.5, 0.0, 0.0], SA[0.0, 0.0, 0.0], 20, :Y, :G),
        ]
    end
    
    @testset "Default values" begin
        d = Dict{String, Any}()
        default_params = Quoll.Parser.PostprocessParams()
        params = from_dict(Quoll.Parser.PostprocessParams, d)

        @test params.fermi_level           == default_params.fermi_level                         
        @test params.dos                   == default_params.dos                                 
        @test params.band_structure        == default_params.band_structure                      
        @test params.fermi_level_params    == default_params.fermi_level_params
        @test params.dos_params            == default_params.dos_params        
        @test params.band_structure_params == default_params.band_structure_params
    end

end