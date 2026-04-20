# Installation

Quoll requires Julia 1.12 or newer.

## As a package

Quoll is not yet registered in the General registry, so it needs to be installed
directly from GitHub. In the Julia REPL, press `]` to enter the package-manager mode
and run:

```julia-repl
pkg> add https://github.com/maurergroup/Quoll.jl.git
```

This gives you access to `using Quoll` from scripts, notebooks, and other packages.

## As an app

If you only need the command-line pipeline, and not the library, install Quoll as a
Julia *app* (a feature introduced in Julia 1.12):

```julia-repl
pkg> app add https://github.com/maurergroup/Quoll.jl.git
```

This makes a `quoll` executable available on your `PATH` (under
`$JULIA_DEPOT_PATH/bin`, which is usually `~/.julia/bin`). You can then run:

```bash
quoll input_file.toml
```

See the [input file documentation](../input_file.md) for the TOML format.

## MPI and HDF5

Quoll is designed to run MPI-parallel pipelines against system MPI and HDF5
libraries — this is almost always what you want on an HPC system. MPI and HDF5
are selected via a `LocalPreferences.toml` file.

When Quoll is installed as an app, the app gets its own environment, typically at:

```
$JULIA_DEPOT_PATH/environments/apps/Quoll/
```

Place a `LocalPreferences.toml` alongside that environment's `Project.toml`.
The repository includes an example
[`LocalPreferences.toml`](https://github.com/maurergroup/Quoll.jl/blob/main/LocalPreferences.toml)
that targets a specific HPC installation of OpenMPI and HDF5:

```toml
[MPIPreferences]
__clear__ = ["preloads_env_switch"]
_format = "1.0"
abi = "OpenMPI"
binary = "system"
cclibs = []
libmpi = "/software/.../OpenMPI/4.1.5-GCC-12.3.0/lib/libmpi.so"
mpiexec = "/software/.../OpenMPI/4.1.5-GCC-12.3.0/bin/mpiexec"
preloads = []

[HDF5]
libhdf5 = "/software/.../HDF5/1.14.0-gompi-2023a/lib/libhdf5.so"
libhdf5_hl = "/software/.../HDF5/1.14.0-gompi-2023a/lib/libhdf5_hl.so"
```

The two sections are independent:

- `[MPIPreferences]` tells [MPI.jl](https://juliaparallel.org/MPI.jl/) which MPI
  implementation to bind to. The same configuration can be produced
  programmatically with `MPIPreferences.use_system_binary()`, which is usually the
  most reliable way to set it up on a new machine.
- `[HDF5]` tells [HDF5.jl](https://juliaio.github.io/HDF5.jl/) which shared
  libraries to use instead of the bundled `HDF5_jll`. The libraries must be ABI
  compatible with the MPI implementation selected above.

!!! note
    After adding `[HDF5]`, `HDF5_jll` may fail to precompile because the prebuilt
    JLL binary no longer matches the libraries you've pointed it at. This is
    expected and harmless — the MPI-linked HDF5 library is still loaded at
    runtime. Precompilation for `Quoll` and `HDF5.jl` themselves should still
    succeed.

Once the preferences are in place, an MPI run looks like:

```bash
mpiexec -n <np> quoll input_file.toml
```

where `<np>` is the number of MPI ranks. Quoll splits work first across input
directories and then across k-points within each directory; see
[What the package does](what_it_does.md) for more detail.
