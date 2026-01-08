module Quoll

using Base
using Reexport
using ArgCheck
using LinearAlgebra
using StaticArrays
using StructArrays
using Memoize
using Dictionaries
using AxisKeys
using Unitful
using DelimitedFiles
using HDF5
using JSON
using AtomsBase
using AtomsIOPython
using NeighbourLists
using Spglib
using MPI

### CONSTANTS ###

include("constants.jl")

### TOOLS ###

include("tools/mpitools.jl")
include("tools/atomtools.jl")
include("tools/fastlookup.jl")
include("tools/symmetry.jl")

### OPERATOR COMPONENTS ###

export Overlap, Hamiltonian
include("core/components/source.jl")
include("core/components/operatorkind.jl")
include("core/components/shconvention.jl")
include("core/components/basis.jl")
include("core/components/spin.jl")
include("core/components/sparsity.jl")

export construct_kgrid
include("core/components/kpoints.jl")

### OPERATOR METHODS AND COMPOSITE TYPES ###

include("core/struct.jl")

export build_operator
include("core/factories.jl")

export find_operatorkinds
export load_operators, load_operator
export write_operators
include("core/inout.jl")

export convert_operator
include("core/conversions.jl")

export fourier_transform
include("core/fourier.jl")

### OPERATORS ###
# Assuming no dependencies between different format types in each file

include("operators/common.jl")
include("operators/canonical.jl")
include("operators/fhiaims.jl")
include("operators/deeph.jl")

### CONVERSIONS ###

include("conversions/canonical.jl")
include("conversions/deeph.jl")

### BASIS PROJECTION MODULE ###

include("Projections.jl")
@reexport using .Projections

### POSTPROCESSING MODULE ###

include("Postprocessing.jl")
@reexport using .Postprocessing

### PARSER MODULE ###

include("Parser.jl")
using .Parser

end
