# Input file

```@meta
CurrentModule = Quoll
```

When Quoll is run as an app (`quoll <input_file.toml>`), its behaviour is driven
entirely by a TOML input file. This page documents every section and key that
is currently wired into the pipeline. An example covering the most common
options is included in the repository at
[`input_file.toml`](https://github.com/maurergroup/Quoll.jl/blob/main/input_file.toml).

!!! note
    A handful of keys are accepted by the parser but are not yet consumed by the
    pipeline (band-structure bands, DOS, Fermi level, error metrics, …). Those
    are intentionally omitted from this page until the corresponding
    functionality lands.

## Top-level structure

```toml
[input]              # required — what to read
[output]             # required — what to write
[kpoint_grid]        # optional — k-point grid used by projection
[basis_projection]   # optional — core-orbital projection
[symmetry]           # optional — symmetry options passed to k-grid construction
```

## `[input]`

Describes the data Quoll should read.

| Key          | Type                 | Default            | Description                                                                         |
| ------------ | -------------------- | ------------------ | ----------------------------------------------------------------------------------- |
| `format`     | string               | —                  | Input format. Currently supported: `"FHI-aims"`, `"FHI-aims-spin"`, `"DeepH"`, `"DeepH-spin"`. |
| `directory`  | string / string list | —                  | One or more root directories. Subdirectories are searched recursively for leaves.   |
| `operators`  | string / string list | all available      | Operators to load. Use `"H"` / `"S"` (or `"Hnospin"`, `"Hspin"`, …) — see below.    |
| `radii`      | string / string list | `nothing`          | Optional per-species interaction radii to build sparsity manually.                  |

### Format names

The format string is case- and whitespace-insensitive. Internally it resolves to
a metadata type via `get_readformat(Val(:...))` — see
[`src/operators/`](https://github.com/maurergroup/Quoll.jl/tree/main/src/operators)
for the mapping.

### Directory

Each directory is used as a root from which Quoll finds every leaf subdirectory
(directories that don't contain further directories). Each leaf is treated as one
atomic configuration. This lets you point Quoll at a tree of many calculations
and have it process all of them in one run.

If multiple root directories are supplied, the output tree for each root is
placed under a subdirectory named after the root's basename, so individual trees
stay separate.

### Operators

The `operators` key accepts shorthand strings that expand to one or more
[`OperatorKind`](@ref)s. The most useful strings are:

| String       | Expands to                                              |
| ------------ | ------------------------------------------------------- |
| `"H"`        | all Hamiltonian variants                                |
| `"S"`        | all Overlap variants                                    |
| `"Hnospin"`  | non-spin Hamiltonians                                   |
| `"Hspin"`    | spin-polarised Hamiltonians (`spin=:up`, `:down`, `:soc`) |
| `"Snospin"`  | non-spin Overlaps                                       |
| `"Sspin"`    | spin-polarised Overlaps                                 |

For bespoke selections, extend the `get_operatorkinds` dispatch table.

### Radii

When provided, each entry has the shape `"<species> <radius> <unit>"`, e.g.
`"Au 5.0 ang"` or `"C 3.0 bohr"`. Two species interact if the distance between
their atoms is at most `r₁ + r₂`. If `radii` is omitted, the sparsity pattern is
taken directly from the input operator data.

## `[output]`

Describes the data Quoll should write.

| Key          | Type     | Default   | Description                                                                  |
| ------------ | -------- | --------- | ---------------------------------------------------------------------------- |
| `format`     | string   | —         | Output format. Currently: `"DeepH"`, `"DeepH-spin"`.                         |
| `directory`  | string   | —         | Output root. The input directory tree is replicated under this root.         |
| `hermitian`  | bool     | inferred  | If `true`, write only the upper triangle; if `false`, write the full matrix. |

If `hermitian` is omitted, Quoll uses the input's own hermicity. Most
DeepH-based training pipelines want `hermitian = false`.

## `[kpoint_grid]`

Controls construction of the k-point grid used by core projection (and, in
future, by postprocessing). Has effect only when an operation that needs a grid
is requested.

| Key        | Type                                    | Default             | Description                                                                                      |
| ---------- | --------------------------------------- | ------------------- | ------------------------------------------------------------------------------------------------ |
| `mesh`     | integer list of length 3                | `nothing`           | Explicit Monkhorst–Pack mesh dimensions. Takes priority over `density`.                          |
| `density`  | float                                   | `10.0`              | Reciprocal-space density in k-points per Å⁻¹, used if `mesh` is not given.                       |
| `shift`    | bool / bool list of length 3            | `[false,false,false]` | Half-grid shift per axis.                                                                        |
| `kpoints`  | string (file path) / list of 4-vectors  | `nothing`           | Explicit list of `[kx, ky, kz, weight]` entries, either inline or read from a whitespace-separated file. Weights should sum to 1. |

Supplying `kpoints` disables `mesh`, `density` and `shift`. When
`kpoints` is a file path, the file is read with `readdlm` — one k-point per line.

For non-periodic systems Quoll always uses a trivial single-Γ-point grid, so
the options above are ignored.

## `[basis_projection]`

Turns on core-orbital projection (see the
[MACE-H / DeepH-E3 tutorial](tutorials/mace_h_deeph_e3.md) for background).
Omit the section to skip projection entirely.

| Key                | Type          | Default   | Description                                                   |
| ------------------ | ------------- | --------- | ------------------------------------------------------------- |
| `projected_basis`  | string list   | —         | Basis subblocks to project out (syntax below).                |
| `method`           | string        | —         | Projection scheme. Currently: `"LaikovCore"`, `"FC99V"`.      |

### Basis syntax

Each string has the shape:

```
<species> (<n>, <l>, <m>) [key = value ...]
```

- `<species>` is the chemical element symbol (e.g. `Au`, `O`, `H`).
- `<n>`, `<l>` are integer principal and angular-momentum quantum numbers.
- `<m>` is either an integer magnetic quantum number or `*` to include all
  `2l + 1` components.
- Any number of `key = value` tags may follow, used to disambiguate multiple
  subblocks that share the same `(n, l, m)`. For FHI-aims data, `type` is the
  most common — the available types are those that appear in
  `basis-indices.out` (e.g. `atomic`, `hydro`, `ionic`).

Example (gold; all core + semicore orbitals):

```toml
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
```

### Requirements

Core projection requires an `Overlap(source=:ref)` among the input operators.

## `[symmetry]`

Tweaks the symmetry reduction used when building the k-point grid.

| Key                | Type       | Default     | Description                                                                 |
| ------------------ | ---------- | ----------- | --------------------------------------------------------------------------- |
| `time_reversal`    | bool       | auto        | Use time-reversal (`k ↔ −k`) to reduce the grid. Auto: on when allowed by the requested operators. |
| `crystal_symmetry` | bool       | auto        | Use crystal point-group symmetry to reduce the grid. Auto: on when there is no projection and the operators allow it. |
| `space_group`      | int        | `nothing`   | Override the detected space group (1–230). Might be used for 2D band structures in the future. |
| `symprec`          | float      | `1e-5`      | Symmetry-detection tolerance in Å passed to Spglib.                         |

Quoll makes conservative automatic choices: spin-polarised operators turn off
the symmetries, and basis-projection turns off crystal symmetry (which the
projection code does not currently handle). You can
force either on by setting the corresponding key explicitly — but see
[`src/parser/methods.jl`](https://github.com/maurergroup/Quoll.jl/blob/main/src/parser/methods.jl)
for the cases where the parser will refuse.
