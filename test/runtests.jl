using Quoll
using Test
using Aqua

@testset "Quoll.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(Quoll)
    end
    # Write your tests here.
end
