using Test

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