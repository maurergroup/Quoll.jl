using SafeTestsets

@safetestset "common.jl" begin
    include("parser/common.jl")
end

@safetestset "input_output.jl" begin
    include("parser/input_output.jl")
end

@safetestset "basis_projection.jl" begin
    include("parser/basis_projection.jl")
end

@safetestset "error_metrics.jl" begin
    include("parser/error_metrics.jl")
end

@safetestset "postprocessing.jl" begin
    include("parser/postprocessing.jl")
end

@safetestset "QuollParams" begin
    include("parser/quollparams.jl")
end