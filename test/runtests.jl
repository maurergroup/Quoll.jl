using Quoll
using Test
using Aqua
using SafeTestsets

@safetestset "Code quality" begin
    include("code_quality.jl")
end

@safetestset "Utils.jl" begin
    include("Utils.jl")
end

@safetestset "Parser.jl" begin
    include("Parser.jl")
end