# Preparing training data for MACE-H / DeepH-E3

This tutorial shows how to use Quoll to convert FHI-aims output into the DeepH
input layout consumed by ML tight-binding training frameworks such as
[MACE-H](https://doi.org/10.1038/s41524-026-02020-1) and
[DeepH-E3](https://github.com/Xiaoxun-Gong/DeepH-E3), with an optional
core-orbital projection step along the way.

It is based on the regression test
[`fhi_canon_laikov_deeph_water.jl`](https://github.com/maurergroup/Quoll.jl/blob/main/test/regression/fhi_canon_laikov_deeph_water.jl),
which exercises the same pipeline on a small water geometry.

## When you'd want a core projection

Deep learning Hamiltonian pipelines typically learn only the *valence* block of
the Hamiltonian — the matrix elements that couple valence basis functions —
because the core states contribute a large but essentially rigid block that is
wasteful to learn. To train on valence-only targets, the core contribution has
to be folded into the valence block rather than simply dropped; otherwise the
eigenvalues of the projected Hamiltonian do not agree with the original valence
eigenvalues.

Quoll implements two such projection schemes, both run by
[`perform_core_projection`](@ref):

- **[`LaikovCore`](@ref)** — the method used in the MACE-H paper. Uses the
  inverse core-core overlap to orthogonalise the core/valence blocks before
  folding.
- **[`FC99V`](@ref)** — a faster frozen-core approximation that avoids the
  inverse; equivalent to `LaikovCore` when the core-core overlap is close to the
  identity.

Both require the reference overlap `S_ref` among the input operators and both
assume hermitian input.

## The input file

The simplest way to run this pipeline is via the app. The TOML below takes
FHI-aims output from `water_fhi/`, projects out the oxygen 1s core orbital, and
writes the valence-only Hamiltonian and overlap in DeepH layout:

```toml
[input]
format = "FHI-aims"
directory = "water_fhi"
operators = ["H", "S"]

[output]
format = "DeepH"
directory = "water_deeph"
hermitian = false

[basis_projection]
projected_basis = [
    "O (1, 0, 0) type = atomic",
]
method = "FC99V"
```

and run it with:

```bash
quoll input_file.toml
```

or, with MPI:

```bash
mpiexec -n <np> quoll input_file.toml
```

The output directory will contain the DeepH files (`hamiltonians.h5`,
`overlaps.h5`, `lat.dat`, `site_positions.dat`, `element.dat`,
`orbital_types.dat`, `R_list.dat`, `rlat.dat`, `info.json`) that MACE-H or
DeepH-E3 read directly.

## The basis syntax

The `projected_basis` field is a list of strings, one per *subblock* of basis
functions to project out. Each entry has the shape:

```
<species> (<n>, <l>, <m>) [key = value ...]
```

where:

- `<species>` is the chemical element symbol.
- `<n>`, `<l>`, `<m>` are the principal, angular-momentum and magnetic quantum
  numbers. Use `*` for `m` to include all `2l+1` magnetic components.
- The optional `key = value` tags distinguish variants of the same `(n, l, m)`
  subblock. For FHI-aims data, `type = atomic` picks the atomic basis as
  reported in `basis-indices.out`; other types (`hydro`, `ionic`, …) can be
  selected in the same way.

A hydrogen atom or a valence-only all-electron basis will usually need no
projection at all, so the entire `[basis_projection]` section can be omitted
for those systems.

## K-point grid

Core projection is performed in reciprocal space (the machinery is the same
whether the system is periodic or not — for an isolated molecule the grid
trivially reduces to the Γ point). The app constructs a k-point grid
automatically. For periodic systems you'll usually want to specify it
explicitly; see the [Input file](../input_file.md) documentation for the full
set of knobs, but a typical section looks like:

```toml
[kpoint_grid]
mesh = [7, 7, 7]
```

## Driving it from the library

If you want more control — e.g. to run just the projection on operators that are
already in memory, or to branch on a custom condition — the same pipeline is
directly accessible from Julia code:

```julia
using Quoll

# Load
H, S = load_operators(
    Operator, FHIaimsCSCNoSpinRealMetadata, "water_fhi",
    [Hamiltonian(source=:ref), Overlap(source=:ref)],
)

# Convert to canonical
H_canon = convert_operator(KeyedOperator, CanonicalBlockRealMetadata, H)
S_canon = convert_operator(KeyedOperator, CanonicalBlockRealMetadata, S)

# Build the k-point grid (single Γ point here, with unit weight)
using StaticArrays
kpoints = [SA[0.0, 0.0, 0.0]]
weights = [1.0]
kgrid = Quoll.KGrid(kpoints, weights, false, false)

# Run the projection
using MPI; MPI.Init()
projected_basis = [
    Quoll.BasisMetadata(
        Quoll.ChemicalSpecies(:O), 1, 0, 0,
        Base.ImmutableDict(:type => :atomic),
    ),
]
H_v, S_v = perform_core_projection(
    [H_canon, S_canon], projected_basis, kgrid, 1:length(kpoints), MPI.COMM_WORLD;
    method = LaikovCore(),
)

# Convert and write
H_out = convert_operator(Operator, DeepHBlockRealMetadata, H_v; hermitian=false)
S_out = convert_operator(Operator, DeepHBlockRealMetadata, S_v; hermitian=false)
write_operators(DeepHBlockRealMetadata, "water_deeph", [H_out, S_out])
```

This is the same flow the app runs — the app is essentially a thin MPI-aware
wrapper around [`perform_core_projection`](@ref), [`convert_operator`](@ref)
and [`write_operators`](@ref).
