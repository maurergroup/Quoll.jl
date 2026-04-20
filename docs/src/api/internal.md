# Internal reference

This page auto-collects docstrings for non-exported ("private") symbols in
`Quoll` and its submodules. They are internal implementation details — useful
to read if you're digging into the code, but not part of the stable API. They
may change without notice between versions.

For the supported public API, see [Public API](public.md); for the extension
interface, see [Operator interface](operator_interface.md).

```@meta
CurrentModule = Quoll
```

## Quoll

```@autodocs
Modules = [Quoll]
Public = false
Private = true
```

## Quoll.Projections

```@autodocs
Modules = [Quoll.Projections]
Public = false
Private = true
```

## Quoll.Postprocessing

```@autodocs
Modules = [Quoll.Postprocessing]
Public = false
Private = true
```

## Quoll.Parser

The `Parser` submodule is not reexported, so all of its symbols are listed here.

```@autodocs
Modules = [Quoll.Parser]
```
