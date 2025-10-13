using Quoll

using Test

@testset "SpinMetadata" begin
    ref = SpinMetadata(1//2)
    @test SpinMetadata(0.5) == ref

    ref = SpinMetadata(-1//2)
    @test SpinMetadata(-0.5) == ref

    @test_throws ArgumentError SpinMetadata(1)
end