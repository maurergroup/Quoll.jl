using SafeTestsets

@safetestset "Utils.jl" begin
    include("Utils.jl")
end

@safetestset "Parser.jl" begin
    include("Parser.jl")
end