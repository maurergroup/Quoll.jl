using Quoll
using Test
using SafeTestsets
using Configurations

include("../testutils.jl")

@safetestset "input_output.jl" begin
    include("Parser/input_output.jl")
end

@safetestset "basis_projection.jl" begin
    include("Parser/basis_projection.jl")
end

@safetestset "error_metrics.jl" begin
    include("Parser/error_metrics.jl")
end

@safetestset "postprocessing.jl" begin
    include("Parser/postprocessing.jl")
end

# TODO: Test all the parameters that were read in once we know the file
# is unlikely to change anymore
@testset "QuollParams" begin
    SiC_exampledir = joinpath(@__DIR__, "..", "..", "examples", "SiC")
    SiC_inputfile = joinpath(SiC_exampledir, "input_file.toml")
    tarballs = [joinpath(SiC_exampledir, "SiC_FHIaims.tar.gz")]

    setupteardown_tmp(tarballs) do
        params = from_toml(Quoll.Parser.QuollParams, SiC_inputfile)
        @test params.input.format == FHIaimsOperator
    end
end