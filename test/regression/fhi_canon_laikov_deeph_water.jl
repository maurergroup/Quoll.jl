# Advanced functionality through the input file, including core orbital projection.
# This regression test is similar to `fhi_canon_laikov_deeph_gold.jl`, but it differs
# because here a non-periodic system is considered, and a different core projection method
# is being used.

using Quoll
using DelimitedFiles
using LazyArtifacts
using StaticArrays
using LinearAlgebra

using Test
using Main.TestUtils

np = 1
project = joinpath(@__DIR__, "../../app")
app = joinpath(@__DIR__, "../../app/quoll.jl")

inputfile = """
[input]
format = "FHI-aims"
directory = "water_fhi/water_fhi"
operators = ["H", "S"]

[output]
format = "DeepH"
directory = "water_deeph_fc99v_converted"
hermitian = false

[basis_projection]
projected_basis = [
    "O (1, 0, 0) type = atomic",
]
method = "FC99V"
"""

test_data = artifact"test_data"
tarballs = [
    joinpath(test_data, "water_fhi.tar.gz")
]

setupteardown_tmp(tarballs = tarballs) do
    inputfile_path = "inputfile.toml"
    write(inputfile_path, inputfile)
    run.(mpiexec_quollapp(app, project, inputfile_path; np=np))

    in_dir = joinpath("water_fhi", "water_fhi")
    out_dir = "water_deeph_fc99v_converted"

    @info "Loading reference eigenvalues"
    eigenvals_ref = readdlm(joinpath(in_dir, "eigenvalues.out"))

    @info "Loading the computed Hamiltonian and overlap operators"
    kind = Quoll.Hamiltonian(source=:ref)
    H = Quoll.load_operator(
        Quoll.Operator, Quoll.DeepHBlockNoSpinRealMetadata, out_dir, kind
    )
    kind = Quoll.Overlap(source=:ref)
    S = Quoll.load_operator(
        Quoll.Operator, Quoll.DeepHBlockNoSpinRealMetadata, out_dir, kind
    )

    @info "Converting operators to canonical format"
    H = Quoll.convert_operator(
        Quoll.KeyedOperator, Quoll.CanonicalBlockRealMetadata, H
    )
    S = Quoll.convert_operator(
        Quoll.KeyedOperator, Quoll.CanonicalBlockRealMetadata, S
    )

    @info "Performing Fourier transform to the gamma point"
    # This is currently the easiest way to generate a dense operator
    kpoint = SVector{3}([0, 0, 0])
    H_dense = pop!(Quoll.fourier_transform(
        Quoll.Operator, Quoll.CanonicalDenseRecipMetadata, H, kpoint
    ))
    S_dense = pop!(Quoll.fourier_transform(
        Quoll.Operator, Quoll.CanonicalDenseRecipMetadata, S, kpoint
    ))

    @info "Computing eigenvalues"
    eigenvals = eigvals(
        Hermitian(Quoll.unwrap_data(Quoll.op_data(H_dense))),
        Hermitian(Quoll.unwrap_data(Quoll.op_data(S_dense)))
    )

    # Compare original eigenvalues to projected, valence-only eigenvalues
    n_core = 2
    Δeigenvals = abs.(eigenvals - eigenvals_ref[begin + n_core : end])
    @test maximum(Δeigenvals) ≈ 0.00013747120267204593
end