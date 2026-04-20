# Public API

This page lists the types and functions exported by `Quoll` and its reexported
submodules. For methods that plug into the extension interface (e.g. when
adding a new format) see [Operator interface](operator_interface.md). For
internal helpers and private functions, see [Internal reference](internal.md).

```@meta
CurrentModule = Quoll
```

## Operators and metadata

```@docs
Operator
KeyedOperator
OperatorKind
```

!!! note
    `Hamiltonian` and `Overlap` are aliases for `OperatorKind{:Hamiltonian}` and
    `OperatorKind{:Overlap}`. Call them with keyword tags: `Hamiltonian(source=:ref)`,
    `Overlap(source=:ref, spin=:up)`, and so on.

## Loading and writing

```@docs
find_operatorkinds
load_operator
load_operators
write_operators
```

## Building and converting

```@docs
build_operator
convert_operator
fourier_transform
```

## K-point grids

```@docs
construct_kgrid
```

## Core-basis projection

```@meta
CurrentModule = Quoll.Projections
```

```@docs
perform_core_projection
LaikovCore
FC99V
```

## Postprocessing

The `Postprocessing` submodule currently exports smearing-function types
(`Gaussian`, `FermiDirac`) used by the input-file parser. They are singleton
types without additional methods attached.

## Auto-indexed

The list below collects any public docstrings not called out individually above.

```@meta
CurrentModule = Quoll
```

```@autodocs
Modules = [Quoll, Quoll.Projections, Quoll.Postprocessing]
Public = true
Private = false
```
