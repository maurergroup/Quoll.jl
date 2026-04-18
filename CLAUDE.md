# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project overview

Quoll.jl is a Julia package for loading, converting, and post-processing tight-binding / DFT operators (Hamiltonians, overlaps) across multiple electronic-structure codes (FHI-aims, DeepH, and a canonical internal format). The package ships a library (`src/`) and a thin MPI-parallel CLI driver (`app/quoll.jl`) that runs a TOML-configured pipeline over directories of operator data.

Julia version: `1.12` (see [Project.toml](Project.toml) compat). The package targets MPI-parallel runs on HPC systems — MPI and HDF5 are typically linked against system libraries via a local `LocalPreferences.toml` (gitignored; must be set up per-machine with `MPIPreferences.use_system_binary` and HDF5 system paths).

## Commands

### Tests

Tests live under [test/](test/) and are gated by CLI flags in [runtests.jl](test/runtests.jl) — running `Pkg.test()` without flags does nothing. Match the CI invocation:

```bash
# Run everything (what CI does):
julia --project -e 'using Pkg; Pkg.test(test_args=["--quality","--unit","--regression"])'

# Subsets:
julia --project -e 'using Pkg; Pkg.test(test_args=["--unit"])'
julia --project -e 'using Pkg; Pkg.test(test_args=["--quality"])'      # Aqua checks
julia --project -e 'using Pkg; Pkg.test(test_args=["--regression"])'   # MPI end-to-end, slow
```

Regression tests download a large test-data artifact (see [Artifacts.toml](Artifacts.toml)) and launch the CLI via `mpiexec` with `np=2` — they require a working MPI install matching the local `LocalPreferences.toml`. Unit tests use `SafeTestsets` and can be drilled into by `include`-ing a single file under [test/unit/](test/unit/) directly in a REPL with `using Main.TestUtils`.

### Running the CLI app

The CLI uses its own environment at [app/](app/) (which `dev`s the parent package). Canonical invocation:

```bash
mpiexec -n <np> julia -t 1 --project=app app/quoll.jl <input_file.toml>
```

TOML input format is driven by `QuollParams` in [src/Parser.jl](src/Parser.jl); the top-level sections are `[input]`, `[output]`, `[kpoint_grid]`, `[basis_projection]`, `[postprocessing]`, `[error_metrics]`, `[symmetry]`. Subparameter schemas live in [src/parser/](src/parser/). The `examples/` directory (gitignored) contains runnable configurations if present on a given machine.

### Formatting

Blue style via [.JuliaFormatter.toml](.JuliaFormatter.toml). Notable non-defaults: `always_use_return = true`, `align_assignment = true`, `align_matrix = true`.

## Architecture

### Pipeline (driven by [app/quoll.jl](app/quoll.jl))

1. Parse TOML → `QuollParams` ([src/Parser.jl](src/Parser.jl)).
2. `find_operatorkinds` + `load_operators` — discover and load operator files in a format-specific layout.
3. `convert_operator` to `CanonicalBlockRealMetadata` — all intra-pipeline work happens in the canonical format.
4. Optional: `construct_kgrid`, `perform_core_projection` (basis projection), postprocessing.
5. Optional: `convert_operator` to an output format and `write_operators`.

Work is split across MPI ranks twice: across input directories (`split_work` with `DefaultSplit`), then across k-points within each config via a sub-communicator (`FHIaimsLAPACKSplit`).

### Operator type system — the core abstraction

Operators are parametrised by a **metadata type** `M <: AbstractMetadata{O,X,S,B,Y,A}` that encodes six things in the type parameters: operator kind, source format, sparsity, basis set, spherical-harmonics convention, and atomic system. See [src/core/struct.jl](src/core/struct.jl). Two Holy traits — `SpaceTrait` (`RealSpace`/`RecipSpace`) and `SpinTrait` (`NoSpin`/`HasSpin`) — dispatch behavior across the four concrete metadata variants (`RealMetadata`, `RecipMetadata`, `SpinRealMetadata`, `SpinRecipMetadata`).

Two operator containers wrap metadata + data: `Operator` and `KeyedOperator` (the latter uses `AxisKeys` for label-based indexing). Conversion between formats is the central operation — always via `convert_operator(OP_type, target_metadata_type, operator)`.

Adding a new format generally means: a file under [src/operators/](src/operators/) defining `load_metadata_basic`, `load_metadata`, `load_data`, `get_avail_operatorkinds`, `get_avail_filenames`, `build_data`, `write_operators`; and a file under [src/conversions/](src/conversions/) implementing `convert_metadata` / `convert_data` between that format's metadata type and the canonical one.

### Module layout

- [src/Quoll.jl](src/Quoll.jl) — top-level module; `include` order matters (constants → tools → components → core struct/factories/inout/conversions/fourier → operator formats → conversions → submodules).
- [src/core/components/](src/core/components/) — building blocks of metadata: `OperatorKind` (kind symbol + tag dict), `BasisMetadata`, `AbstractSparsity`, `AbstractSource`, `SHConvention`, spin, k-points.
- [src/operators/](src/operators/) — format-specific loaders/writers/data builders (`canonical`, `fhiaims`, `deeph`).
- [src/conversions/](src/conversions/) — metadata/data converters to/from the canonical format.
- [src/tools/](src/tools/) — `mpitools` (`split_work`, `DefaultSplit`, `FHIaimsLAPACKSplit`), `atomtools`, `fastlookup`, `symmetry`.
- Submodules reexported through `Quoll`:
  - `Projections` ([src/Projections.jl](src/Projections.jl)) — `FC99V`, `LaikovCore` core-state basis projectors; entry point `perform_core_projection`.
  - `Postprocessing` ([src/Postprocessing.jl](src/Postprocessing.jl)) — smearings (`Gaussian`, `FermiDirac`).
  - `Parser` (not reexported) — TOML input parsing, assembled into `QuollParams` with validation hooks (`validate_operatorkinds`, `validate_symmetry`, `search_clashes`).

### Units and conventions

Internal computations use `UnitfulAtomic`. Lengths at I/O boundaries use `LengthAngstrom` (see [src/constants.jl](src/constants.jl)). Magnetic quantum numbers `m` are always stored in the "wiki" SH convention internally; `SHConvention` tracks the source/target convention for conversions.
