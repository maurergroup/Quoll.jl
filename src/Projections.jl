module Projections

using LinearAlgebra
using StaticArrays
using ArgCheck
using MPI

using ..Quoll
using ..Quoll:
    AbstractOperator,
    KGrid,
    precompute_phases,
    get_float,
    get_kind,
    get_data,
    get_metadata,
    get_sparsity,
    get_basisset,
    get_spins,
    get_dense_subbasis_mask,
    construct_subbasis_metadata,
    build_operator,
    fourier_transform,
    inv_fourier_transform_data!,
    synchronise_data!

abstract type AbstractBasisProjection end

function get_basis_projection end

function perform_core_projection end

export perform_core_projection
include("projections/common.jl")

export FC99V
include("projections/fc99v.jl")

export LaikovCore
include("projections/laikov.jl")

end