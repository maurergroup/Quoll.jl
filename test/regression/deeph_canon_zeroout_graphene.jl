# Regression test for `zero_out_data!`.
#
# Two Hamiltonians describe the *same* graphene system, holding (up to numerical noise) the
# same matrix elements, but with different sparsity:
#
#   * `graphene_deeph`        — a MACE-H prediction in DeepH format, stored on a *denser*
#                               neighbour list.
#   * `graphene_fhi_injected` — the same prediction injected into FHI-aims and dumped right
#                               after injection, i.e. the data projected onto FHI-aims' own
#                               (sparser) CSC sparsity pattern.
#
# `zero_out_data!` restricts an operator to a sparser support without changing its type: it
# zeroes every entry of the (denser) DeepH operator whose sparsity is not contained in the
# FHI-aims CSC sparsity, while leaving the shared entries untouched. Because the FHI-aims
# operator is by construction the DeepH data projected onto that same CSC sparsity, the two
# operators must agree once the DeepH operator has been zeroed out.
#
# The test therefore checks:
#   1. before zeroing, the two canonical operators differ (DeepH carries a few extra
#      long-range sub-blocks that FHI-aims truncated);
#   2. after zeroing the DeepH operator against the FHI-aims CSC metadata, they agree.

using Quoll
using StaticArrays
using Dictionaries
using LazyArtifacts

using Test
using Main.TestUtils

# Tolerance for comparing the two Hamiltonians. The shared (non-zeroed) matrix elements are
# only equal up to the numerical error introduced by FHI-aims' internal FT/invFT loop during
# injection (see the README of the data). Tune if the regression flags a spurious mismatch.
compare_kwargs = (atol = 1e-8, rtol = 1e-6)

# Flatten a canonical `KeyedOperator` into a Dict keyed by `(iat, jat, R)` → dense block, so
# that two operators with *different* sparsity can be compared element-wise: a block absent
# from an operator's sparsity is treated as an implicit zero block by `same_data` below.
function collect_blocks(operator)
    keydata = Quoll.unwrap_data(Quoll.op_keydata(operator))
    ij2images = Quoll.op_sparsity(operator).ij2images
    blocks = Dict{Tuple{Int,Int,SVector{3,Int}},Matrix{Float64}}()
    for ((iat, jat), block) in pairs(keydata)
        images = ij2images[(iat, jat)]
        for k in eachindex(images)
            blocks[(iat, jat, images[k])] = Matrix{Float64}(collect(block[:, :, k]))
        end
    end
    return blocks
end

# Do two canonical operators carry the same data, treating any block missing from either
# operator's sparsity as a zero block?
function same_data(op1, op2; kwargs...)
    blocks1 = collect_blocks(op1)
    blocks2 = collect_blocks(op2)
    for key in union(keys(blocks1), keys(blocks2))
        block1 = get(blocks1, key, nothing)
        block2 = get(blocks2, key, nothing)
        block1 = isnothing(block1) ? zero(block2) : block1
        block2 = isnothing(block2) ? zero(block1) : block2
        isapprox(block1, block2; kwargs...) || return false
    end
    return true
end

test_data = artifact"test_data"
tarballs = [
    joinpath(test_data, "graphene_deeph.tar.gz")
    joinpath(test_data, "graphene_fhi_injected.tar.gz")
]

setupteardown_tmp(tarballs = tarballs) do
    deeph_dir = joinpath("graphene_deeph", "graphene_deeph")
    fhi_dir = joinpath("graphene_fhi_injected", "graphene_fhi_injected")

    @info "Loading the DeepH (denser) and FHI-aims-injected (sparser) Hamiltonians"
    H_deeph = load_operator(
        Quoll.Operator, Quoll.DeepHBlockNoSpinRealMetadata, deeph_dir,
        Quoll.Hamiltonian(source=:pred),
    )
    H_fhi = load_operator(
        Quoll.Operator, Quoll.FHIaimsCSCNoSpinRealMetadata, fhi_dir,
        Quoll.Hamiltonian(source=:ref),
    )
    @test Quoll.op_metadata(H_deeph) isa Quoll.DeepHBlockNoSpinRealMetadata
    @test Quoll.op_metadata(H_fhi) isa Quoll.FHIaimsCSCNoSpinRealMetadata

    # Keep the FHI-aims CSC metadata: its (sparser) sparsity is what we zero the DeepH
    # operator against, and it must be captured before `H_fhi` is converted to canonical.
    fhi_csc_metadata = Quoll.op_metadata(H_fhi)

    @info "Converting both operators to the canonical block format"
    # The DeepH prediction is dumped with full (non-hermitian) storage, so hermiticity is
    # forced on conversion — `zero_out_data!` only supports the hermitian-to-hermitian case.
    H_deeph_canon = convert_operator(
        Quoll.KeyedOperator, Quoll.CanonicalBlockRealMetadata, H_deeph; hermitian=true
    )
    H_fhi_canon = convert_operator(
        Quoll.KeyedOperator, Quoll.CanonicalBlockRealMetadata, H_fhi
    )
    @test Quoll.op_metadata(H_deeph_canon) isa Quoll.CanonicalBlockNoSpinRealMetadata
    @test Quoll.op_metadata(H_fhi_canon) isa Quoll.CanonicalBlockNoSpinRealMetadata

    @info "Before zeroing out: the operators should differ (DeepH carries extra blocks)"
    @test !same_data(H_deeph_canon, H_fhi_canon; compare_kwargs...)

    @info "Zeroing out the DeepH operator against the FHI-aims CSC sparsity"
    Quoll.zero_out_data!(H_deeph_canon, fhi_csc_metadata)

    @info "After zeroing out: the operators should agree"
    @test same_data(H_deeph_canon, H_fhi_canon; compare_kwargs...)
end
