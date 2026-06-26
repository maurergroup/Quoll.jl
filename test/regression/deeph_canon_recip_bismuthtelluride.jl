# Regression test covering a spin (SOC) operator and a full Fourier-transform round-trip.
#
# The Hamiltonian is forward-transformed onto a dense k-point grid and inverse-transformed
# back to real space, then converted to DeepH and written out. The written `hamiltonians.h5`
# is compared against the original input, which must be recovered to floating-point precision.
#
# This exercises the complex-output branch of `inv_fourier_transform_data!`: the SOC operator
# has `ComplexF64` data, so the inverse transform must keep the complex contribution (a full,
# unreduced grid is used). Dropping the imaginary part would break the round-trip.

using Quoll
using StaticArrays
using LazyArtifacts
using HDF5

using Test
using Main.TestUtils

test_data = artifact"test_data"
tarballs = [
    joinpath(test_data, "bismuthtelluride_deeph.tar.gz")
]

setupteardown_tmp(tarballs = tarballs) do
    dir = joinpath("bismuthtelluride_deeph", "bismuthtelluride_deeph")
    operatorkind = Quoll.Hamiltonian(source=:ref, spin=:soc)

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

    @info "Constructing a dense k-point grid (no symmetry reduction for a spin operator)"
    atoms = Quoll.op_atoms(H_canonical)
    kgrid = construct_kgrid(
        atoms; density=10.0, time_reversal=false, crystal_symmetry=false
    )
    symmetry = Quoll.grid_symmetry(kgrid)
    @test symmetry == Quoll.KGridSymmetry{false,false}()

    kpoints_float = convert(Vector{SVector{3,Float64}}, kgrid.kpoints)
    phases = collect(
        eachcol(Quoll.precompute_phases(kpoints_float, Quoll.op_images(Quoll.op_sparsity(H_canonical))))
    )

    @info "Fourier transform → inverse Fourier transform round-trip over the grid"
    H_recovered = build_operator(
        Quoll.KeyedOperator, Quoll.op_metadata(H_canonical);
        value=zero(Quoll.op_float_type(H_canonical)), initialised=true,
    )

    for (ik, (kpoint, weight)) in enumerate(zip(kpoints_float, kgrid.weights))
        H_recip = fourier_transform(
            Quoll.KeyedOperator, Quoll.CanonicalDenseRecipMetadata, H_canonical, kpoint,
            phases[ik]; out_shconv=Quoll.op_shconv(H_canonical),
        )
        Quoll.inv_fourier_transform_data!(
            H_recovered, H_recip, phases[ik], weight, symmetry
        )
    end

    @info "Converting recovered real-space operator back to DeepH and writing"
    H_out = convert_operator(Quoll.Operator, Quoll.DeepHBlockSpinRealMetadata, H_recovered)
    out_dir = "bismuthtelluride_deeph_roundtrip"
    write_operators(Quoll.DeepHBlockSpinRealMetadata, out_dir, [H_out])

    @testset "hamiltonians.h5" begin
        db_ref = h5open(joinpath(dir, "hamiltonians.h5"))
        db = h5open(joinpath(out_dir, "hamiltonians.h5"))

        k_ref = sort(keys(db_ref))
        k = sort(keys(db))

        blocks_ref = read.(Ref(db_ref), k_ref)
        blocks = read.(Ref(db), k)

        close(db_ref)
        close(db)

        @test k == k_ref
        @test blocks ≈ blocks_ref
    end
end
