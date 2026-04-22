# Quoll

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://maurergroup.github.io/Quoll.jl/dev/)
[![Build Status](https://github.com/maurergroup/Quoll.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/maurergroup/Quoll.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Aqua](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

Quoll.jl is a Julia package for loading, converting, and post-processing
quantum-mechanics operators (Hamiltonians, overlaps, ...) expressed in an
atomic-orbital basis. It provides a single internal *canonical* representation
that acts as a pivot between code-specific formats such as
[FHI-aims](https://fhi-aims.org/),
[DeepH](https://github.com/mzjb/DeepH-pack), and
[MACE-H](https://github.com/maurergroup/MACE-H),
and it handles the bookkeeping around such operators — sparsity patterns, basis
sets, spherical-harmonics conventions, k-point grids, and spin channels.

A typical use case is producing training data for machine-learning electronic-structure models:
for example, taking real-space Hamiltonians and overlaps from
an FHI-aims calculation, projecting out the core basis functions, and writing
the result in the DeepH layout expected by
[DeepH-E3](https://github.com/Xiaoxun-Gong/DeepH-E3) or
[MACE-H](https://doi.org/10.1038/s41524-026-02020-1)
training pipelines.

## Features

- **Multiple formats.** Built-in support for FHI-aims, DeepH, and Quoll's own
  canonical format. Conversion between any pair is a single `convert_operator`
  call.
- **Extensible via dispatch.** Adding a new format, sparsity pattern, or
  operator flavour does not require forking Quoll — define the relevant methods
  on your own metadata type, and the rest of the pipeline works unchanged.
- **Core-basis projection** onto minimal core bases (`FC99V`, `LaikovCore`).
- **Fourier transforms** between real- and reciprocal-space representations on
  user-specified k-point grids.
- **MPI-parallel CLI app** that runs a TOML-configured pipeline over
  directories of operator data and over k-points on HPC systems.

## Installation

Quoll requires Julia 1.12 or newer and is not yet registered in the General
registry.

As a library:

```julia-repl
pkg> add https://github.com/maurergroup/Quoll.jl.git
```

As a command-line app (Julia 1.12+):

```julia-repl
pkg> app add https://github.com/maurergroup/Quoll.jl.git
```

This installs a `quoll` executable under `$JULIA_DEPOT_PATH/bin` (usually
`~/.julia/bin`).

On HPC systems, you will typically want to bind MPI.jl and HDF5.jl to the
system libraries via a `LocalPreferences.toml`; see the
[installation docs](https://maurergroup.github.io/Quoll.jl/dev/getting_started/installation/)
for details.

## Quick start

### As a library

```julia
using Quoll

dir   = "examples/water/water_fhi"
kinds = find_operatorkinds(Quoll.FHIaimsCSCNoSpinRealMetadata, dir)
H, S  = load_operators(Quoll.Operator, Quoll.FHIaimsCSCNoSpinRealMetadata, dir, kinds)

H_canon = convert_operator(Quoll.KeyedOperator, Quoll.CanonicalBlockRealMetadata, H)
S_canon = convert_operator(Quoll.KeyedOperator, Quoll.CanonicalBlockRealMetadata, S)

H_deeph = convert_operator(Quoll.Operator, Quoll.DeepHBlockRealMetadata, H_canon; hermitian=false)
S_deeph = convert_operator(Quoll.Operator, Quoll.DeepHBlockRealMetadata, S_canon; hermitian=false)

write_operators(Quoll.DeepHBlockRealMetadata, "water_deeph", [H_deeph, S_deeph])
```

### As a CLI app

Describe the same pipeline in a TOML file:

```toml
[input]
format    = "FHI-aims"
directory = "water_fhi"
operators = ["H", "S"]

[output]
format    = "DeepH"
directory = "water_deeph"
hermitian = false
```

and run:

```bash
quoll input_file.toml
# or, in parallel over MPI ranks:
mpiexec -n <np> quoll input_file.toml
```

Ready-to-run examples live under [examples/](examples/). For more complicated cases, take a look at the regression tests.

## Documentation

Full documentation — installation notes, the TOML input-file reference,
tutorials (including MACE-H / DeepH-E3 training-data preparation), the public
API, and guidance on extending Quoll — is available at

<https://maurergroup.github.io/Quoll.jl/dev/>

## License

Quoll.jl is released under the [MIT License](LICENSE).
