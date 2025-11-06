using Quoll
using Configurations

using Test
using Main.TestUtils

@testset "normalize_comparison" begin

    s = "FERMI_DIRAC"
    @test Quoll.Parser.normalize_comparison(s) == "fermidirac"

    replace_pairs = tuple()
    @test Quoll.Parser.normalize_comparison(s; replace_pairs) == "fermi_dirac"

end