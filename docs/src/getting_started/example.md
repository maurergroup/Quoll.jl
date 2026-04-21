# First example

```@meta
CurrentModule = Quoll
DocTestSetup = quote
    using Quoll
    cd(joinpath(pkgdir(Quoll), "examples", "water"))
    rm("water_deeph"; force=true, recursive=true)
end
DocTestTeardown = quote
    rm(joinpath(pkgdir(Quoll), "examples", "water", "water_deeph");
       force=true, recursive=true)
    cd(joinpath(pkgdir(Quoll), "docs"))
end
```

This example walks through the smallest end-to-end workflow in Quoll: read a
Hamiltonian and overlap from an FHI-aims calculation, convert them to the
canonical format, then convert to DeepH and write to disk. Along the way we
introduce the main abstractions you need to know about.

The data comes from `examples/water/water_fhi/` in the Quoll repository — an
FHI-aims calculation for two isolated water molecules in a non-periodic box.
All commands below assume you've `cd`'d into `examples/water/`.

The full list of names is in the [API reference](../api/core.md); everything
non-exported is qualified with `Quoll.` throughout this page.

## The library version

```julia
using Quoll
```

### Loading

Quoll represents every loaded matrix as an **operator** with a tagged
**metadata** type. The metadata type encodes six pieces of information through
its type parameters: the operator *kind* (Hamiltonian vs. overlap), the source
format, the sparsity pattern, the basis set, the spherical-harmonics convention,
and the atomic structure.

The first thing we need to do is discover which operators FHI-aims actually
produced in our directory:

```jldoctest pipeline
julia> dir = "water_fhi";

julia> kinds = find_operatorkinds(Quoll.FHIaimsCSCNoSpinRealMetadata, dir)
2-element Vector{Quoll.OperatorKind}:
 Hamiltonian({:source = :ref})
 Overlap({:source = :ref})
```

[`find_operatorkinds`](@ref) scans `dir` for any files that a given metadata
type knows how to read (here `rs_hamiltonian.out`, `rs_overlap.out`, etc.) and
returns the list of `OperatorKind`s it found. An `OperatorKind` is a
little more than a symbol: it's a kind (`:Hamiltonian` or `:Overlap`) plus a
dictionary of tags like `source=:ref` or `spin=:up`, which disambiguate variants
that would otherwise have the same name. `Hamiltonian` and `Overlap` above are
exported aliases for `OperatorKind{:Hamiltonian}` and
`OperatorKind{:Overlap}`.

Now we can load them:

```jldoctest pipeline
julia> H, S = load_operators(Quoll.Operator, Quoll.FHIaimsCSCNoSpinRealMetadata, dir, kinds);

julia> H
Operator(RealMetadata(Hamiltonian({:source = :ref}), FHIaimsSource, CSCRealSparsity))
```

[`load_operators`](@ref) dispatches on the metadata type
(`FHIaimsCSCNoSpinRealMetadata`) to pick the right file-reading code, and
returns concrete `Operator`s — the simplest of the two operator
containers. `KeyedOperator` is the richer alternative; it additionally
wraps the matrix data in per-atom-pair views indexed by `(i, j)` pairs, which
is what most manipulation code in Quoll uses internally.

### Conversion to the canonical format

All intra-pipeline work happens in the canonical format. We convert with
[`convert_operator`](@ref):

```jldoctest pipeline
julia> H_canon = convert_operator(Quoll.KeyedOperator, Quoll.CanonicalBlockRealMetadata, H)
KeyedOperator(RealMetadata(Hamiltonian({:source = :ref}), CanonicalSource, BlockRealSparsity))

julia> S_canon = convert_operator(Quoll.KeyedOperator, Quoll.CanonicalBlockRealMetadata, S);
```

`convert_operator` does three things under the hood:

1. **Metadata conversion** (`convert_metadata`): change the source, the
   spherical-harmonics convention, the sparsity pattern, and (optionally) reduce
   the basis set.
2. **Allocation** ([`build_operator`](@ref)): allocate an empty output operator
   with zero-initialised data.
3. **Data conversion** (`convert_data!`): transfer and transform the
   actual matrix elements.

The first argument selects the output operator *container* (`Operator` or
`KeyedOperator`), the second selects the output *metadata* type. Note
that the metadata type here is a **union alias** —
`CanonicalBlockRealMetadata` covers both the spin-polarised and the
non-spin variants; the concrete subtype is chosen automatically based on
whether the input has spin information.

### Conversion to DeepH and writing

The conversion to DeepH follows the same pattern:

```jldoctest pipeline
julia> H_deeph = convert_operator(Quoll.Operator, Quoll.DeepHBlockRealMetadata, H_canon; hermitian=false)
Operator(RealMetadata(Hamiltonian({:source = :ref}), DeepHSource, BlockRealSparsity))

julia> S_deeph = convert_operator(Quoll.Operator, Quoll.DeepHBlockRealMetadata, S_canon; hermitian=false);
```

The `hermitian=false` keyword asks Quoll to write the full matrix rather than
only the upper triangle — this is what DeepH-based training pipelines typically
expect.

Finally, write the result:

```jldoctest pipeline
julia> write_operators(Quoll.DeepHBlockRealMetadata, "water_deeph", [H_deeph, S_deeph])

julia> sort(readdir("water_deeph"))
9-element Vector{String}:
 "R_list.dat"
 "element.dat"
 "hamiltonians.h5"
 "info.json"
 "lat.dat"
 "orbital_types.dat"
 "overlaps.h5"
 "rlat.dat"
 "site_positions.dat"
```

[`write_operators`](@ref) dispatches on the metadata type to pick the right
file layout. For DeepH this produces `hamiltonians.h5`, `overlaps.h5`, and the
associated `lat.dat`, `site_positions.dat`, `element.dat`, `orbital_types.dat`,
`R_list.dat`, `rlat.dat` and `info.json` files.

## The app version

The same workflow expressed as a TOML file (`input_file.toml`):

```toml
[input]
format = "FHI-aims"
directory = "water_fhi"
operators = ["H", "S"]

[output]
format = "DeepH"
directory = "water_deeph"
hermitian = false
```

Run it with:

```bash
quoll input_file.toml
```

(no `mpiexec` for this example — only a single configuration is being converted in this case,
and no operations are using k-points.)

This file is what ships in `examples/water/input_file.toml`; it's the shortest
input file the app will accept. The full list of sections and keys is documented
in [Input file](../input_file.md).

## Where to go from here

- To see how core-basis projection, k-point grids, and MPI fit in, read the
  [MACE-H / DeepH-E3 tutorial](../tutorials/mace_h_deeph_e3.md).
- To learn how to plug in a new format or sparsity pattern without modifying
  Quoll.jl, read [Defining new methods](defining_new_methods.md).

```@meta
DocTestSetup = nothing
DocTestTeardown = nothing
```
