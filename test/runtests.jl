using SafeTestsets

# Include TestUtils module here; it can then be used
# across separate test files via `using Main.TestUtils`
include("TestUtils.jl")

@safetestset "Code quality" begin
    include("code_quality.jl")
end

@safetestset "Unit tests" begin
    include("unit/unittests.jl")
end

@safetestset "Regression tests" begin
    include("regression/regressiontests.jl")
end
