using SafeTestsets

@safetestset "common.jl" begin
    include("parser/common.jl")
end

@safetestset "basis_projection.jl" begin
    include("parser/basis_projection.jl")
end

@safetestset "error_metrics.jl" begin
    include("parser/error_metrics.jl")
end

@safetestset "input_output.jl" begin
    include("parser/input_output.jl")
end

@safetestset "kpoint_grid.jl" begin
    include("parser/kpoint_grid.jl")
end

@safetestset "postprocessing.jl" begin
    include("parser/postprocessing.jl")
end

@safetestset "symmetry.jl" begin
    include("parser/symmetry.jl")
end

@safetestset "methods.jl" begin
    include("parser/methods.jl")
end