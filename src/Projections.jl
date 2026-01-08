module Projections

using LinearAlgebra
using StaticArrays
using ArgCheck
using MPI

using ..Quoll
using ..Quoll:
    AbstractOperator,
    AbstractMetadata,
    Operator,
    CanonicalDenseNoSpinRecipMetadata,
    KeyedTrait,
    SpinTrait,
    KGrid,
    DenseRecipData,
    op_float_type,
    op_data_type,
    op_kind,
    op_metadata,
    op_data,
    op_sparsity,
    op_basisset,
    op_shconv,
    op_hermicity,
    op_images,
    precompute_phases,
    get_subbasis_masks,
    get_dense_subbasis_mask,
    wrap_data,
    unwrap_data,
    convert_metadata,
    build_operator,
    build_operator_extra,
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