# Regression test that covers operators with spin, and an application of Fourier transform.
# Currently this test does not compare the values obtained, but it could be extended.

using Quoll
using StaticArrays
using LazyArtifacts

using Test
using Main.TestUtils

test_data = artifact"test_data"
tarballs = [
    joinpath(test_data, "bismuthtelluride_deeph.tar.gz")
]

setupteardown_tmp(tarballs = tarballs) do
    dir = joinpath("bismuthtelluride_deeph", "bismuthtelluride_deeph")
    operatorkind = Quoll.Hamiltonian(source=:ref, spin=:soc)
    kpoint = SVector{3}([0, 0.5, 0.5])

    @info "Loading DeepH Operator"
    H_deeph = load_operator(
        Quoll.Operator, Quoll.DeepHBlockSpinRealMetadata, dir, operatorkind
    )
    @test Quoll.op_metadata(H_deeph) isa Quoll.DeepHBlockSpinRealMetadata

    @info "Converting to Canonical KeyedOperator"
    H_canonical = convert_operator(
        Quoll.KeyedOperator, Quoll.CanonicalBlockRealMetadata, H_deeph
    )
    @test Quoll.op_metadata(H_canonical) isa Quoll.CanonicalBlockSpinRealMetadata

    @info "Performing Fourier transform to compute reciprocal-space dense KeyedOperator"
    H_canonical_recip = pop!(fourier_transform(
        Quoll.KeyedOperator, Quoll.CanonicalDenseRecipMetadata, H_canonical, kpoint
    ))
    @test Quoll.op_metadata(H_canonical_recip) isa Quoll.CanonicalDenseSpinRecipMetadata
end
