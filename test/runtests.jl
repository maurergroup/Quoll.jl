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

# TODO: This could be changed to Pkg.test() julia ARG instead
const JULIA_QUOLL_TEST_REGRESSION = parse(Bool, get(ENV, "JULIA_QUOLL_TEST_REGRESSION", "true"))

if JULIA_QUOLL_TEST_REGRESSION
    @safetestset "Regression tests" begin
        include("regression/regressiontests.jl")
    end
end