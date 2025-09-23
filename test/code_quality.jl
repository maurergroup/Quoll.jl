using Quoll
using Test
using Aqua

@testset "Code Quality (Aqua.jl)" begin
    Aqua.test_ambiguities(Quoll)
    Aqua.test_unbound_args(Quoll)
    Aqua.test_undefined_exports(Quoll)
    Aqua.test_stale_deps(Quoll)
    Aqua.test_deps_compat(
        Quoll,
        check_extras = false,
        check_weakdeps = false)
    Aqua.test_piracies(Quoll)
    Aqua.test_undocumented_names(Quoll)
end
