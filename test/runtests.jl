using SafeTestsets

@safetestset "Code quality" begin
    include("code_quality.jl")
end

@safetestset "Unit tests" begin
    include("unit/unittests.jl")
end

@safetestset "Regression tests" begin
    include("regression/regressiontests.jl")
end