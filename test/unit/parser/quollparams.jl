using Quoll
using Configurations

using Test
using Main.TestUtils

# TODO: Test all the parameters that were read in once we know the file
# is unlikely to change anymore
@testset "SiC" begin
    SiC_exampledir = joinpath(@__DIR__, "../../../examples/SiC")
    SiC_inputfile = joinpath(SiC_exampledir, "input_file.toml")
    tarballs = [joinpath(SiC_exampledir, "SiC_FHIaims.tar.gz")]

    setupteardown_tmp(tarballs = tarballs) do
        params = from_toml(Quoll.Parser.QuollParams, SiC_inputfile)
        @test params.input.format == FHIaimsCSCOperator
    end
end