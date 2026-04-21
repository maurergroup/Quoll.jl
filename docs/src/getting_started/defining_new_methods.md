# Defining new methods

```@meta
CurrentModule = Quoll
```

Quoll is built around multiple dispatch on **metadata types**: a metadata type
encodes the operator kind, source format, sparsity pattern, basis set, SH
convention and atomic structure all in its type parameters, and every
format-specific method in the pipeline keys off one of those parameters. As a
result most extensions — new formats, new sparsity patterns, new operator
flavours, new conversion paths — can be written **outside** of Quoll.jl.
As long as the methods are on the method table before the pipeline runs,
multiple dispatch picks them up.

This page walks through what a typical extension looks like. For the precise
method signatures, see [Operator interface](../api/operator_interface.md); for
the existing format implementations used as references throughout, see
[`src/operators/fhiaims.jl`](https://github.com/maurergroup/Quoll.jl/blob/main/src/operators/fhiaims.jl),
[`src/operators/deeph.jl`](https://github.com/maurergroup/Quoll.jl/blob/main/src/operators/deeph.jl),
and
[`src/operators/canonical.jl`](https://github.com/maurergroup/Quoll.jl/blob/main/src/operators/canonical.jl).

## Defining a new format

A "format" in Quoll is the combination of a source tag, a metadata alias,
and a data layout.

### Source singleton

Every format is identified by a subtype of [`AbstractSource`](@ref):

```julia
struct MyFormatSource <: Quoll.AbstractSource end
```

If your format uses a non-standard spherical-harmonics convention (sign
conventions, ordering of `m` components), define its `default_shconv`:

```julia
@memoize function Quoll.default_shconv(::MyFormatSource)
    return Quoll.SHConvention(...)
end
```

### Data alias

Operator data is wrapped in a [`DataContainer`](@ref) parametrised by element
type, dimensionality, body container type, source, and sparsity. You declare
your format's layout by writing a `const` alias, the same way
`CSCRealData` and `DenseRecipData` are declared in
[`src/core/struct.jl`](https://github.com/maurergroup/Quoll.jl/blob/main/src/core/struct.jl):

```julia
const MyFormatData{T} = Quoll.DataContainer{
    T,                         # floating point type
    N,                         # dimensionality
    <:SomeBodyType{T},         # underlying array or dictionary
    <:MyFormatSource,
    <:SomeSparsity,
}
```

The body type here is whatever makes sense for your storage — a single
`Array{T,N}`, a per-atom-pair `Dictionary`, etc. The source and sparsity types
are included as well — this allows to further differentiate data formats in
multiple dispatch when `SomeBodyType` is the same, but has different structure,
e.g. two body types could be dictionaries, but have different dictionary key
conventions not captured by the type.

### Metadata aliases

A metadata alias specialises one of the concrete metadata variants
([`RealMetadata`](@ref), [`RecipMetadata`](@ref),
[`SpinRealMetadata`](@ref), [`SpinRecipMetadata`](@ref)) by pinning its
source and sparsity type parameters. Writing aliases for each (spin, space)
combination you support mirrors what
`FHIaimsCSCNoSpinRealMetadata`, `FHIaimsCSCSpinRealMetadata` and the
union alias `FHIaimsCSCRealMetadata` do in
[`src/operators/fhiaims.jl`](https://github.com/maurergroup/Quoll.jl/blob/main/src/operators/fhiaims.jl):

```julia
const MyFormatNoSpinRealMetadata{
    O<:OperatorKind,
    X<:MyFormatSource,
    S<:SomeSparsity,
    B<:BasisSetMetadata,
    Y<:SHConvention,
    A<:AbstractSystem
} = Quoll.RealMetadata{O,X,S,B,Y,A}

const MyFormatSpinRealMetadata{
    O<:OperatorKind,
    X<:MyFormatSource,
    S<:SomeSparsity,
    B<:BasisSetMetadata,
    Y<:SHConvention,
    A<:AbstractSystem,
    P<:SpinsMetadata
} = Quoll.SpinRealMetadata{O,X,S,B,Y,A,P}

const MyFormatRealMetadata{...} = Union{
    <:MyFormatNoSpinRealMetadata{...},
    <:MyFormatSpinRealMetadata{...},
}
```

The union alias is what most of the machinery dispatches on; the concrete
no-spin / spin aliases are used where the distinction actually matters (e.g.
`load_metadata`).

### Type-level hooks

Three single-line methods tell the pipeline how to query types from a
metadata type:

```julia
Quoll.op_data_type(::Type{<:MyFormatRealMetadata})     = MyFormatData
Quoll.op_source_type(::Type{<:MyFormatRealMetadata})   = MyFormatSource
Quoll.op_sparsity_type(::Type{<:MyFormatRealMetadata}) = SomeSparsity
```

If you want [`KeyedOperator`](@ref) support, also declare a keydata alias
and point to it with `op_keydata_type`:

```julia
Quoll.op_keydata_type(::Type{<:MyFormatRealMetadata}) = MyFormatKeyData
```

### Defining a new metadata type

The four concrete metadata variants cover every (real/recip space) ×
(no-spin/has-spin) combination. You only need a new concrete metadata
**type** if your format carries a field none of them currently do.

When that happens, the new type needs:

- a `trait(::Type{T}, ::Type{<:MyMetadata})` method for each trait it
  participates in (e.g. [`SpaceTrait`](@ref), [`SpinTrait`](@ref));
- an `extrafield_traittypes` method listing which traits control
  its extra fields, so that [`convert_metadata`](@ref)'s extras pass and
  [`build_operator`](@ref)'s extras pass know what to build;
- conversion handlers for the new trait — see the next section.

For a worked example, [`src/core/struct.jl`](https://github.com/maurergroup/Quoll.jl/blob/main/src/core/struct.jl)
shows how `SpaceTrait` and `SpinTrait` are attached to the four existing
metadata types.

## Two ways to get metadata into existence

Once the metadata aliases are in place, the pipeline produces concrete
metadata values along one of two paths: loading from disk, or converting from
another metadata value. You implement whichever matches how the user will use
your format.

### The `load_metadata` path

[`load_operator`](@ref) calls `load_metadata` followed by `load_data`,
both dispatched on your metadata alias. The default `load_metadata` is a
two-stage function: it calls `load_metadata_basic` to assemble a
`BasicMetadataContainer`, then calls a second `load_metadata` method that
wraps it into the correct concrete metadata type. You define both:

```julia
function Quoll.load_metadata_basic(
    ::Type{<:MyFormatRealMetadata}, dir, kind::Quoll.OperatorKind,
)
    source   = MyFormatSource()
    atoms    = load_my_atoms(dir)
    sparsity = load_my_sparsity(dir, kind)
    basisset = load_my_basisset(dir, atoms)
    shconv   = Quoll.default_shconv(source)
    return Quoll.BasicMetadataContainer(kind, source, sparsity, basisset, shconv, atoms)
end

function Quoll.load_metadata(
    ::Type{<:MyFormatNoSpinRealMetadata}, dir, basic::Quoll.BasicMetadataContainer,
)
    return Quoll.RealMetadata(basic)
end
```

You can optionally also declare which `OperatorKind`s your format supports on disk and
which files back each one (used by the app):

```julia
Quoll.get_avail_operatorkinds(::Type{<:MyFormatNoSpinRealMetadata}) = [
    Hamiltonian(; source=:ref),
    Overlap(; source=:ref),
]

Quoll.get_avail_filename(
    ::Type{<:MyFormatNoSpinRealMetadata},
    ::Tuple{Val{:Hamiltonian},Pair{Val{:source},Val{:ref}}},
) = "my_hamiltonian.dat"
```

These are what [`find_operatorkinds`](@ref) scans with — the default
implementation checks, for each `OperatorKind` returned by
`get_avail_operatorkinds`, whether any of its `get_avail_filenames`
exists on disk.

### The `convert_metadata` path

If someone wants to convert another format's metadata into your format,
[`convert_operator`](@ref) will call [`convert_metadata`](@ref), which calls:

1. [`convert_metadata_basic`](@ref) — converts source, SH convention,
   sparsity and basis set. This stage is generic and needs no per-format
   method unless your sparsity type is also new (in which case you implement
   a format-specific [`convert_sparsity`](@ref) for it).
2. `convert_metadata_extra` — dispatches on each trait in
   `extrafield_traittypes(Mₒᵤₜ)` to convert extra fields (k-point, spins).
   The generic methods for `SpaceTrait` and `SpinTrait` cover the four
   existing metadata types already; you only write new
   `convert_metadata_extra` methods if you introduced a new trait.
3. `convert_metadata_final` — assembles the concrete output metadata type
   from the basic container and the extras. Again, only needed for a new
   metadata type.

For the common case of "new format, same metadata variants," step 1 is all
that happens automatically and neither of the `_extra`/`_final` hooks should
need touching.

## Building the operator

[`build_operator`](@ref) takes a metadata value and allocates a matching
operator. The default implementation handles `Operator` and `KeyedOperator`
generically — the two things you typically need to provide are data and,
for keyed operators, keydata.

### `build_data`

[`build_data`](@ref) is dispatched on `op_data_type(M)`, so one method
per data alias:

```julia
function Quoll.build_data(
    ::Type{<:MyFormatData}, metadata::M, value::T, initialised::Bool,
) where {M<:Quoll.AbstractMetadata,T<:Number}
    body = allocate_my_body(metadata, T)
    initialised && fill_body!(body, value)
    return Quoll.wrap_data(M, body)
end
```

`wrap_data` tags the raw container with the source and sparsity types so
that dispatch elsewhere in the pipeline knows what it is.

### `build_keydata`

Only needed if you want [`KeyedOperator`](@ref) support. Given the
just-allocated data, produce the keyed view of it:

```julia
function Quoll.build_keydata(
    ::Type{<:MyFormatKeyData}, metadata, data::MyFormatData,
)
    return Quoll.wrap_data(typeof(metadata), build_my_keyed_view(metadata, data))
end
```

### `build_operator_extra` / `build_operator_final`

You only need these if you introduce a **new operator type** (a new
subtype of `AbstractOperator` that carries extra fields of its own),
and they mirror the `convert_metadata_extra` / `convert_metadata_final`
pattern: each trait in `extrafield_traittypes(OP)` gets a hook that
builds its extra field, and `build_operator_final` assembles the
concrete operator. For `Operator` / `KeyedOperator` the existing
implementations cover everything you need.

## Converting the operator

Data conversion is the point at which actual matrix elements are moved
between formats. Metadata conversion has already happened by the time the
data-level conversion methods are invoked, so your method receives fully
constructed input and output operators and is responsible only for the
transfer.

### `convert_data!`

Defined per data-type pair. Which signature to implement depends on the
[`KeyedTrait`](@ref) of each operator:

```julia
# Neither operator is keyed
function Quoll.convert_data!(
    ::Type{<:MyFormatData}, ::Type{<:SomeOtherData},
    out_op::Quoll.Operator, in_op::Quoll.Operator,
)
    # transfer in_op.data to out_op.data
end

# Output is keyed, input isn't (and the three other combinations)
function Quoll.convert_data!(
    ::Type{<:MyFormatKeyData}, ::Type{<:MyFormatData}, ::Type{<:SomeOtherData},
    out_op, in_op,
)
    ...
end
```

Adding a new format typically means one method per direction, between your
data type and the canonical data type — `src/conversions/` in Quoll.jl is
structured exactly that way.

### `fourier_transform_data!` and `inv_fourier_transform_data!`

These follow the same 4-way `KeyedTrait` dispatch as `convert_data!`,
with an extra `phases_k` argument (and, for the inverse transform, a
`weight` for k-point summation). Implement them only if your format
participates in Fourier transforms — typically this means a real-space
format converting to a dense reciprocal-space format, or vice versa. If
your format is real-space only (or reciprocal-space only) and the Fourier
path runs through the canonical format, you can skip these.

## Operator I/O

Loading and writing are optional — you can define a format that exists
purely as an in-memory intermediate, or as a conversion target with no
on-disk representation. They are only needed if the user will feed a
directory of your format into the pipeline or ask Quoll to write it.

### Loading

[`load_operator`](@ref) dispatches to your `load_metadata` +
`load_metadata_basic` (covered above) and then `load_data`:

```julia
function Quoll.load_data(
    ::Type{M}, dir, kind::Quoll.OperatorKind,
) where {M<:MyFormatRealMetadata}
    body = read_my_body(dir, kind)
    return Quoll.wrap_data(M, body)
end
```

The file names `load_data` opens should match what
`get_avail_filenames` / `get_avail_filename` report; this
is what keeps [`find_operatorkinds`](@ref), [`load_operator`](@ref) and
[`write_operators`](@ref) in agreement with each other.

### Writing

[`write_operators`](@ref) is declared as an abstract `function ... end`
stub in [`src/core/inout.jl`](https://github.com/maurergroup/Quoll.jl/blob/main/src/core/inout.jl); each format provides its own concrete method:

```julia
function Quoll.write_operators(
    ::Type{M}, dir, operators,
) where {M<:MyFormatRealMetadata}
    ispath(dir) || mkpath(dir)
    for op in operators
        write_my_data(dir, op)
    end
    write_my_shared_metadata(dir, first(operators))
end
```

Same consistency rule as above: the files you produce should match the
names advertised by `get_avail_filenames`, so that the same directory
round-trips back through `find_operatorkinds` + `load_operator`.

## Driving it from the TOML app

Once the library-level methods exist, two one-line registrations expose
the format to the CLI:

```julia
Quoll.get_readformat(::Val{:myformat})  = MyFormatNoSpinRealMetadata
Quoll.get_writeformat(::Val{:myformat}) = MyFormatNoSpinRealMetadata
```

After that, `format = "myformat"` in an `[input]` or `[output]` section of
an `input_file.toml` is resolved to your metadata type by
[`src/parser/input_output.jl`](https://github.com/maurergroup/Quoll.jl/blob/main/src/parser/input_output.jl).

## When to upstream

If the extension is specific to your project (a one-off file format, a
custom sparsity trick), keep it in your own code. If it's a generally
useful format or conversion — a new electronic-structure code, a widely
used machine-learning data layout, a standard core-projection scheme —
please consider contributing it. The [developer docs](../developer.md)
describe how to set up a development checkout and what is expected of a
contribution.
