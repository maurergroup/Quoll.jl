# Defining new methods

Quoll is designed so that most extensions — new formats, new sparsity patterns,
new operator flavours, new conversion paths — can be written **outside** of
Quoll.jl and still plug into the rest of the pipeline via multiple dispatch.
This page is an overview of what those extensions look like; see
[Operator interface](../api/operator_interface.md) for the precise method
signatures.

## Is everything optional?

Almost. Which methods you need depends entirely on what you want to do:

| If you want to...                                           | You need...                                                              |
| ----------------------------------------------------------- | ------------------------------------------------------------------------ |
| Load a new format from disk                                 | `load_metadata_basic`, `load_metadata`, `load_data`, `get_avail_filenames`, `get_avail_operatorkinds` |
| Write a new format to disk                                  | `write_operators`                                                        |
| Convert from/to your format                                 | `convert_metadata` and `convert_data!` (for every format pair you care about) |
| Allocate zero-initialised operators of your format          | `build_data`, `op_data_type`, `op_source_type`, `op_sparsity_type`       |
| Use `KeyedOperator` with your format                        | `build_keydata`, `op_keydata_type`                                       |
| Drive it from the TOML app                                  | a `get_readformat(::Val{:myformat})` (and `get_writeformat` for output)  |
| Project out a custom core-basis subset                      | nothing extra — the existing `perform_core_projection` works unchanged   |
| Add a new core-projection scheme (like `LaikovCore`)        | an `AbstractBasisProjection` subtype plus a `compute_valence_data` method, and `get_basis_projection(::Val{:myscheme})` if you want it TOML-driven |

So it's genuinely possible, for example, to define just `convert_metadata` and
`convert_data!` between your format and the canonical format if you already
have loaded operators in memory and only need a conversion path. Or to define
only `load_metadata_basic`, `load_metadata` and `load_data` if you're
extending an existing format with a new file layout.

## The shape of a typical extension

Suppose you want to add support for a new file format `MyFormat`. In your own
package (or scratch file), you would:

1. Define a source singleton:
   ```julia
   struct MyFormatSource <: Quoll.AbstractSource end
   ```
2. Define metadata *aliases* for each (spin, space) combination you support —
   analogous to how `FHIaimsCSCNoSpinRealMetadata` is defined as a specialised
   `RealMetadata` in [`src/operators/fhiaims.jl`](https://github.com/maurergroup/Quoll.jl/blob/main/src/operators/fhiaims.jl).
3. Give that alias the data-type, source-type and sparsity-type metadata hooks:
   ```julia
   Quoll.op_data_type(::Type{<:MyFormatMetadata})     = MyFormatData
   Quoll.op_source_type(::Type{<:MyFormatMetadata})   = MyFormatSource
   Quoll.op_sparsity_type(::Type{<:MyFormatMetadata}) = ...
   ```
4. Implement the methods you actually need — loading, writing, conversion —
   from the table above.

The key observation is that **none of this needs to live inside Quoll.jl**. As
long as the methods are on the method table before you call the pipeline,
multiple dispatch will pick them up.

## When to upstream

If the extension is specific to your project (a one-off file format, a custom
sparsity trick), keep it in your own code. If it's a generally useful format or
conversion — a new electronic-structure code, a widely used machine-learning
data layout, a standard core-projection scheme — please consider contributing
it. The [developer docs](../developer.md) describe how to set up a development
checkout and what is expected of a contribution.
