```@meta
CurrentModule = Quoll
```

# Quoll.jl

Quoll is a Julia package for loading, converting and post-processing tight-binding
and DFT operators (Hamiltonians, overlaps, ...) that are expressed in an atomic-orbital
basis. It provides a single internal *canonical* representation that acts as a pivot
between code-specific formats such as
[FHI-aims](https://fhi-aims.org/) and
[DeepH](https://github.com/mzjb/DeepH-pack),
and it handles the bookkeeping that surrounds such operators — sparsity patterns,
basis sets, spherical-harmonics conventions, k-point grids and spin channels.

Quoll ships in two forms:

- a **library** (`using Quoll`) that exposes the type system and conversion pipeline for
  use from scripts, notebooks, or other packages, and
- a **CLI app** (`quoll <input_file.toml>`) that runs a TOML-configured pipeline
  over directories of operator data, optionally in parallel via MPI.

A typical use case is producing training data for machine-learning tight-binding
models — for example, taking real-space Hamiltonians and overlaps from an FHI-aims
calculation, projecting out the core basis functions, and writing the result in the
DeepH layout expected by
[DeepH-E3](https://github.com/Xiaoxun-Gong/DeepH-E3) or
[MACE-H](https://doi.org/10.1038/s41524-026-02020-1)
training pipelines. A full walkthrough of that workflow is given in the
[tutorial](tutorials/mace_h_deeph_e3.md).

## Where to start

- New users should read the [getting-started guide](getting_started/installation.md),
  which covers installation and walks through a small water example.
- For the TOML interface used by the app, see [Input file](input_file.md).
- For the types and functions exposed by the library, see the
  [Public API](api/public.md) and [Operator interface](api/operator_interface.md)
  pages.
- Contributors should look at the [developer docs](developer.md).
