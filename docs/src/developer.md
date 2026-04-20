# Developer docs

This page is for people who want to contribute to Quoll.jl itself. For writing
downstream extensions that don't modify Quoll's source tree, see
[Defining new methods](getting_started/defining_new_methods.md) and
[Operator interface](api/operator_interface.md) instead.

## Development checkout

Clone the repository and install it — plus the app — in development mode:

```bash
git clone https://github.com/maurergroup/Quoll.jl.git
cd Quoll.jl
```

```julia-repl
pkg> dev .
pkg> app add .
```

`dev .` makes the package editable: changes to `src/` are picked up on the next
`using Quoll` without reinstalling. `app add .` additionally registers a
development copy of the `quoll` executable; this is only needed if you plan to
run the regression tests locally, since they shell out to `mpiexec quoll`.

If you are on a machine with system MPI / HDF5, also create the
`LocalPreferences.toml` described in [Installation](getting_started/installation.md).

## Contributing

Contributions are made via pull requests against `main`. A typical PR:

1. Adds or changes code under `src/`.
2. Adds matching tests under `test/unit/` or `test/regression/`.
3. Runs cleanly under the formatter (`.JuliaFormatter.toml`).

Nothing here is rigidly enforced, but tests in particular are appreciated — any
new dispatch method or conversion path is likely to need at least a unit test
to prevent future regressions.

### Extending vs. contributing

If the new functionality is format- or project-specific, there's no need to
upstream it — Quoll is designed so that new methods can be defined outside of
Quoll.jl and still plug into the pipeline. See
[Defining new methods](getting_started/defining_new_methods.md).

If it is of general interest — a new electronic-structure code, a standard ML
data layout, a commonly useful projection scheme — please open an issue or PR.

## Tests

Tests are organised into three groups, gated by CLI flags in
[`test/runtests.jl`](https://github.com/maurergroup/Quoll.jl/blob/main/test/runtests.jl).
`Pkg.test()` with no flags will do nothing, so always pass the groups you
actually want to run:

```bash
# What CI runs:
julia --project -e 'using Pkg; Pkg.test(test_args=["--quality","--unit","--regression"])'

# Unit tests only (fast):
julia --project -e 'using Pkg; Pkg.test(test_args=["--unit"])'

# Aqua-style code quality checks:
julia --project -e 'using Pkg; Pkg.test(test_args=["--quality"])'

# End-to-end regression tests (slow, needs MPI):
julia --project -e 'using Pkg; Pkg.test(test_args=["--regression"])'
```

### Adding unit tests

Unit tests live under [`test/unit/`](https://github.com/maurergroup/Quoll.jl/tree/main/test/unit)
and are organised to mirror the `src/` layout. Each file is a `@safetestset`
included from `test/unit/unittests.jl`. To add coverage:

1. Put the new file under the corresponding subdirectory
   (e.g. `test/unit/core/convert_metadata.jl` for tests of
   `src/core/conversions.jl`).
2. Add an `@safetestset` entry for it in `test/unit/unittests.jl`.
3. Use the helpers from `test/TestUtils.jl` (`using Main.TestUtils`) for
   temporary-directory setup, MPI-free logging, and so on. These can also be
   used to include a single test file directly in a REPL for iteration —
   `using Main.TestUtils; include("test/unit/.../foo.jl")`.

### Adding regression tests

Regression tests live under
[`test/regression/`](https://github.com/maurergroup/Quoll.jl/tree/main/test/regression)
and run the `quoll` CLI end-to-end via `mpiexec`. Each one:

1. Writes out a TOML input file to a temporary directory.
2. Extracts reference data from the `test_data` artifact (see
   [`Artifacts.toml`](https://github.com/maurergroup/Quoll.jl/blob/main/Artifacts.toml)).
3. Runs the app with `mpiexec_quollapp(...)` from `TestUtils.jl`.
4. Compares the produced operators against a reference (e.g. eigenvalues from
   FHI-aims).

New regression tests follow the same pattern. Tooling for updating the
`test_data` artifact (uploading a new tarball, refreshing its hash in
`Artifacts.toml`) is not fully automated yet, so for now just flag in the PR
that you have new reference data and the maintainers will help wire it in.

## Formatting

Quoll uses [JuliaFormatter](https://domluna.github.io/JuliaFormatter.jl/) with
a Blue-style base and a few overrides in
[`.JuliaFormatter.toml`](https://github.com/maurergroup/Quoll.jl/blob/main/.JuliaFormatter.toml)
(notably `always_use_return`, `align_assignment`, `align_matrix`). You can run
it locally with:

```julia-repl
pkg> activate --temp
pkg> add JuliaFormatter
julia> using JuliaFormatter; format(".")
```

Fine-grained formatting details aren't strictly enforced, but sticking to the
configured defaults keeps diffs small.
