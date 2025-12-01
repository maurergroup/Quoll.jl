module Quoll

using Base
using Reexport
using ArgCheck
using LinearAlgebra
using StaticArrays
using OffsetArrays
using Dictionaries
using AutoHashEquals
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

# Constants

include("constants.jl")

# Methods

export recentre
include("methods/mpitools.jl")
include("methods/atomtools.jl")
include("methods/fastlookup.jl")
include("methods/symmetry.jl")

# Core types and their methods

export Overlap, Hamiltonian
include("operatorkind.jl")

export BasisMetadata
include("basis.jl")

export ⬆, ⬇
include("spin.jl")
include("sparsity.jl")
include("shconversion.jl")

export get_kgrid
include("kpoints.jl")

# Operator formats and their IO routines
# (Assuming no dependencies between different format types in each file)

export load_atoms
export load_operator, load_operators
export write_operator, write_operators
export find_operatorkinds
include("operators/interface.jl")

export BSparseOperator, DenseOperator
include("operators/canonical/abstract.jl")
include("operators/canonical/bsparse.jl")
include("operators/canonical/dense.jl")

export DeepHOperator
include("operators/deeph/deeph.jl")

export FHIaimsCSCOperator
include("operators/fhiaims/abstract.jl")
include("operators/fhiaims/fhiaims_csc.jl")

# Operator format conversions

export convert_operator, convert_operators
export fourier_transform
include("conversions/interface.jl")

include("conversions/canonical/bsparse.jl")
include("conversions/fhiaims/fhiaims_csc.jl")
include("conversions/deeph/deeph.jl")

# Basis projection module

include("Projections.jl")
@reexport using .Projections

# Postprocessing module

include("Postprocessing.jl")
@reexport using .Postprocessing

# Parser module

include("Parser.jl")
using .Parser

end
