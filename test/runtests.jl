using SafeTestsets

# Include TestUtils module here; it can then be used
# across separate test files via `using Main.TestUtils`
include("TestUtils.jl")

if "--quality" in ARGS
    @safetestset "Code quality" begin
        include("code_quality.jl")
    end
end

if "--unit" in ARGS
    @safetestset "Unit tests" begin
        include("unit/unittests.jl")
    end
end

if "--regression" in ARGS
    @safetestset "Regression tests" begin
        include("regression/regressiontests.jl")
    end
end