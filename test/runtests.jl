using Quoll
using Test
using Aqua
using SafeTestsets

@safetestset "Code Quality" begin
    include("code_quality.jl")
end
