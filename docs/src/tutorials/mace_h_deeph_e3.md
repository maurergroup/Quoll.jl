# Preparing training data for MACE-H / DeepH-E3

This tutorial shows how to use Quoll to convert FHI-aims output into the DeepH
input layout consumed by ML Hamiltonian training frameworks such as
[MACE-H](https://github.com/maurergroup/MACE-H) and
[DeepH-E3](https://github.com/Xiaoxun-Gong/DeepH-E3), with a core-orbital
projection step along the way.

It is based on the regression test
[`fhi_canon_laikov_deeph_gold.jl`](https://github.com/maurergroup/Quoll.jl/blob/main/test/regression/fhi_canon_laikov_deeph_gold.jl),
which exercises the same pipeline on a bulk-gold geometry — the generated data
matches the data used in the original [MACE-H paper](https://doi.org/10.1038/s41524-026-02020-1).

## When you'd want a core projection

If only valence electronic properties are of interest, then projecting out the
core contribution from the Hamiltonian might be beneficial, as
the core states contribute a large but essentially rigid block that is
wasteful to learn. To train on valence-only targets, the core contribution has
to be folded into the valence block rather than simply dropped; otherwise the
eigenvalues of the projected Hamiltonian do not agree with the original valence
eigenvalues.

Quoll implements two such projection schemes, both run by
[`perform_core_projection`](@ref Quoll.Projections.perform_core_projection):

- **[`LaikovCore`](@ref)** — the method used in the MACE-H paper. Uses the
  inverse core-core overlap to orthogonalise the core/valence blocks before
  folding.
- **[`FC99V`](@ref)** — a faster frozen-core approximation that avoids the
  inverse; equivalent to `LaikovCore` when the core-core overlap is close to the
  identity.

Both require the reference overlap `S_ref` among the input operators and both
assume hermitian input.

## The input file

The TOML below takes FHI-aims output from `gold_fhi/` (can be obtained from test-data branch),
builds sparsity from a 5 Å per-species interaction radius, projects out all Au core and semicore
orbitals with `LaikovCore`, and writes the valence-only Hamiltonian and overlap in DeepH layout:

```toml
[input]
format = "FHI-aims"
directory = "gold_fhi"
operators = ["H", "S"]
radii = ["Au 5.0 ang"]

[output]
format = "DeepH"
directory = "gold_deeph"
hermitian = false

[kpoint_grid]
mesh = [7, 7, 7]

[basis_projection]
projected_basis = [
    "Au (1, 0, 0) type = atomic",
    "Au (2, 0, 0) type = atomic",
    "Au (3, 0, 0) type = atomic",
    "Au (4, 0, 0) type = atomic",
    "Au (5, 0, 0) type = atomic",
    "Au (2, 1, *) type = atomic",
    "Au (3, 1, *) type = atomic",
    "Au (4, 1, *) type = atomic",
    "Au (5, 1, *) type = atomic",
    "Au (3, 2, *) type = atomic",
    "Au (4, 2, *) type = atomic",
    "Au (4, 3, *) type = atomic",
]
method = "LaikovCore"
```

Run it with:

```bash
mpiexec -n <np> quoll input_file.toml
```

Work is split across MPI ranks by k-point, so for the `[7, 7, 7]` mesh above
anywhere up to the (symmetry-reduced) k-point count scales well. For a
local test (might take around 10 minutes as it involves taking the inverse of
the core-core overlap matrix):

```bash
mpiexec -n 4 quoll input_file.toml
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

The twelve subblocks above are every Au orbital up to and including the 5p
semicore — exactly the set frozen out in the MACE-H paper. For a lighter
element, or a valence-only all-electron basis, the list is correspondingly
shorter; if no projection is needed at all, omit the `[basis_projection]`
section entirely.

## Interaction radii

`radii` is how you control the sparsity pattern of the output operators. Each
entry has the shape `"<species> <radius> <unit>"`; two species are considered
to interact when the distance between their atoms is at most `r₁ + r₂`. For
the MACE-H gold example, a single 5 Å radius per Au atom was found to yield
negligible errors compared to the full radius ($\leq$6 Å).

If `radii` is omitted, Quoll takes the sparsity directly from the input
operator data.

## K-point grid

Core projection is performed in reciprocal space, so a k-point grid is
required. `[kpoint_grid]` is where you control it; a Gamma-centered grid
is the default choice:

```toml
[kpoint_grid]
mesh = [7, 7, 7]
```

See [Input file](../input_file.md) for the full set of knobs (explicit
k-points, density instead of mesh, half-grid shift, etc.).

For an isolated molecule the grid trivially reduces to the Γ point regardless
of the settings above — the periodic machinery still works, it just degenerates.
