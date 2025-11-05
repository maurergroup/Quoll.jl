module Quoll

using Reexport
using LinearAlgebra
using StaticArrays
using Dictionaries
using AutoHashEquals
using AxisKeys
using Unitful
using HDF5
using AtomsBase
using NeighbourLists
using ArgCheck

# Common tools

export recentre
include("common/constants.jl")
include("common/mpitools.jl")
include("common/atomtools.jl")
include("common/fastlookup.jl")

# Core types and their methods

export Overlap, Hamiltonian
include("operatorkind.jl")

export BasisMetadata
include("basis.jl")

export ⬆, ⬇
include("spin.jl")
include("sparsity.jl")
include("shconversion.jl")

# Operator formats and their IO routines
# (Assuming no dependencies between different format types in each file)

export load_atoms
export load_operator, load_operators
export write_operator, write_operators
export find_operatorkinds
include("operators/interface.jl")

export RealBSparseOperator
include("operators/canonical/abstract.jl")
include("operators/canonical/bsparse.jl")

export DeepHOperator
include("operators/deeph/deeph.jl")

export FHIaimsCSCOperator
include("operators/fhiaims/abstract.jl")
include("operators/fhiaims/fhiaims_csc.jl")

# Operator format conversions

export convert_operator, convert_operators
include("conversions/interface.jl")

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
