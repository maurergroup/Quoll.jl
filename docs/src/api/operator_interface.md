# Operator interface

This page lists the dispatch points that a new format, sparsity pattern, or
conversion path can hook into. The surrounding narrative for *how* to extend
Quoll is in
[Defining new methods](../getting_started/defining_new_methods.md); this page
is the reference.

All methods below live in the `Quoll` namespace. They are not `export`ed, but
they are stable and may be specialised by downstream code without touching
Quoll.jl itself.

```@meta
CurrentModule = Quoll
```

## Metadata hooks (required per format)

A format is identified by its metadata type (e.g. `FHIaimsCSCRealMetadata`).
Every format must resolve these type-level queries:

| Function              | Returns                                                  |
| --------------------- | -------------------------------------------------------- |
| `op_data_type(::Type{M})`     | The `DataContainer` alias for the operator's data |
| `op_source_type(::Type{M})`   | The `AbstractSource` singleton type               |
| `op_sparsity_type(::Type{M})` | The `AbstractSparsity` subtype                    |
| `op_keydata_type(::Type{M})`  | The `DataContainer` alias for keyed data (needed only if your format works with `KeyedOperator`) |

See [`src/operators/fhiaims.jl`](https://github.com/maurergroup/Quoll.jl/blob/main/src/operators/fhiaims.jl)
for a minimal example.

## Loading

To load operators from disk in your format, implement:

| Function                                        | Purpose                                                  |
| ----------------------------------------------- | -------------------------------------------------------- |
| `load_metadata_basic(::Type{M}, dir, kind)`     | Build a `BasicMetadataContainer` from files on disk      |
| `load_metadata(::Type{M}, dir, basic_metadata)` | Wrap the basic container into the concrete metadata type (Real/Recip/SpinReal/SpinRecip) |
| `load_data(::Type{M}, dir, kind)`               | Read the data array from disk                            |
| `get_avail_filenames(::Type{M}, kind)`          | Return candidate on-disk filenames for this operator kind |
| `get_avail_operatorkinds(::Type{M})`            | Return the list of operator kinds this format supports  |

`load_operator` (the entry point used by the pipeline) is already implemented
generically — it composes the methods above.

```@docs; canonical=false
load_operator
load_operators
find_operatorkinds
```

## Writing

Implement `write_operators(::Type{M}, dir, operators)` for your metadata type.
For multi-file layouts it is usually convenient to split the implementation
into metadata-writing and data-writing helpers — see
[`src/operators/deeph.jl`](https://github.com/maurergroup/Quoll.jl/blob/main/src/operators/deeph.jl)
for the shape.

```@docs; canonical=false
write_operators
```

## Building

For [`build_operator`](@ref) to produce zero-initialised operators in your format:

```@docs; canonical=false
build_data
```

`build_keydata(::Type{KD}, metadata, data)` is additionally required if your
format should work with `KeyedOperator`.

## Converting

To support conversion between your format and another, dispatch on the pair of
metadata types (input and output):

```@docs; canonical=false
convert_metadata
convert_metadata_basic
convert_data!
convert_sparsity
```

In practice, most conversions between your format and the rest of Quoll will
go through `CanonicalBlockRealMetadata`, so the pair to implement is usually
`(YourMetadata, CanonicalBlockRealMetadata)` and vice versa. Direct conversions
to any other format are also supported — just add the method for that pair.

## TOML parsing

If you want your format to be selectable from `input_file.toml`, register it
with the parser's format lookup table:

```julia
Quoll.get_readformat(::Val{:myformat}) = MyFormatMetadata
Quoll.get_writeformat(::Val{:myformat}) = MyFormatMetadata
```

The string after `Val(` is matched against the `format` field after
case/whitespace/punctuation normalisation, so `"MyFormat"`, `"my-format"`, and
`"my_format"` all resolve to the same symbol.

Similarly, custom core-projection schemes register themselves via:

```julia
Quoll.Projections.get_basis_projection(::Val{:myscheme}) = MyScheme
```

## Holy traits

Dispatch on behavioural traits rather than concrete types is used in a handful
of places (real vs. reciprocal space, spin vs. no-spin, keyed vs. unkeyed
operator).

```@docs; canonical=false
SpaceTrait
SpinTrait
KeyedTrait
```
