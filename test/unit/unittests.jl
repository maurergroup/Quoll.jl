using SafeTestsets

@safetestset "Utils.jl" begin
    include("Utils.jl")
end

@safetestset "Parser.jl" begin
    include("Parser.jl")
end

@safetestset "MPITools.jl" begin
    include("MPITools.jl")
end

@safetestset "AtomsTools" begin
    include("AtomsTools.jl")
end