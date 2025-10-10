module Quoll

using Reexport
using LinearAlgebra
using StaticArrays
using Dictionaries
using AxisKeys
using Unitful
using HDF5
using AtomsBase
using ArgCheck

# Common tools
export recentre
include("common/mpitools.jl")
include("common/atomtools.jl")

# Core types and their methods
export BasisMetadata
export Overlap, Hamiltonian
include("basis.jl")
include("operatorkind.jl")
include("sparsity.jl")

# Operator formats and their IO routines
# (Assuming no dependencies between different format types in each file)
export load_atoms, load_operators, find_operatorkinds
include("formats/abstract.jl")

include("formats/canonical/abstract.jl")
include("formats/canonical/bsparse.jl")

include("formats/deeph/deeph.jl")

export FHIaimsCSCOperator
include("formats/fhiaims/abstract.jl")
include("formats/fhiaims/fhiaims_csc.jl")

# Operator format conversions
include("conversions/canonical/bsparse.jl")
include("conversions/fhiaims/fhiaims_csc.jl")

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
