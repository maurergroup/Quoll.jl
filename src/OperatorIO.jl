# TODO: should probably define an API to make it easier for developers
# to implement new formats 
module OperatorIO
using StaticArrays
using AtomsBase
using Dictionaries
using AxisKeys
using ArgCheck
using HDF5
using DelimitedFiles

export read_format, write_format, operator_kind
export avail_operatorkinds
export find_operatorkinds
export load_operators, load_operator_data, load_operator_metadata

export AbstractOperator
export FHIaimsOperator, DeepHOperator
export FHIaimsOperatorMetadata
export AbstractOperatorKind
export Hamiltonian, Overlap

abstract type AbstractOperator end

abstract type AbstractOperatorKind end

struct Hamiltonian <:AbstractOperatorKind
    tag::Symbol
end

struct Overlap <:AbstractOperatorKind
    tag::Symbol
end

abstract type AbstractOperatorMetadata end

abstract type AbstractQuollOperator <: AbstractOperator end
abstract type AbstractBSparseOperator <: AbstractQuollOperator end

# TODO: could the same be used for density matrices from FHI-aims?
# If not, and if I needed to load in H, S, D at the same time it could
# become tricky
struct FHIaimsOperatorMetadata <: AbstractOperatorMetadata
    n_basis::Int
    row_ind::Vector{Int}
    col_cell_ptr::Array{Int, 3}
    cells::Vector{SVector{3, Int}}
    # spherical harmonics convention?
end

# TODO: Same as for metadata, if we want to read in density matrices or gradients,
# we would likely need different fields but the format is still "FHIaims".
# If it's stored not as CSC but something else we could e.g. use parametric
# types instead of using Vector{T}
struct FHIaimsOperator{O<:AbstractOperatorKind, T<:AbstractFloat} <: AbstractOperator
    kind::O
    data::Vector{T}
    metadata::FHIaimsOperatorMetadata
end

struct DeepHOperator <: AbstractOperator end

struct BSparseOperatorMetadata <: AbstractOperatorMetadata
    z1z2_ij2interval::Dictionary{NTuple{2, ChemicalSpecies}, Dictionary{NTuple{2, Int64}, UnitRange{Int64}}}
    atom2species::Vector{ChemicalSpecies}
    # TODO: I think I can keep those three as separate structs and only keep
    # attributes which are required for matrix operations in this struct
    # Sparsity?
    # Atoms?
    # BasisSetMetadata?
    # TODO: possibly construct this by only taking BasisSetMetadata, Atoms, Sparsity,
    # but create an outer constructor which could compute z1z2_ij2interval from those internally
    # TODO: also possibly create a constructor which could create a matrix with zero values
end

struct RealBSparseOperator{O<:AbstractOperatorKind, T<:AbstractFloat, AT, KT} <: AbstractBSparseOperator
    kind::O
    data::Dictionary{NTuple{2, ChemicalSpecies}, Array{T, 3}}
    keydata::Dictionary{NTuple{2, Int}, KeyedArray{T, 3, AT, KT}}
    metadata::BSparseOperatorMetadata
end

function read_format end
function write_format end

read_format(::Val{:fhiaims}) = FHIaimsOperator
write_format(::Val{:deeph}) = DeepHOperator

operator_kind(::Val{:h}) = Hamiltonian(:ref)
operator_kind(::Val{:href}) = Hamiltonian(:ref)
operator_kind(::Val{:hpred}) = Hamiltonian(:pred)

operator_kind(::Val{:s}) = Overlap(:ref)
operator_kind(::Val{:sref}) = Overlap(:ref)
operator_kind(::Val{:spred}) = Overlap(:pred)

function load_operator_metadata(dir::AbstractString, ::Type{FHIaimsOperator})
    @debug "Loading operator metadata"

    p = joinpath(dir, "rs_indices.out")
    @argcheck ispath(p)
    
    f = open(p, "r")

    # Subtracting 1 to ignore one element (last element)
    # which is always zero
    line = readline(f)
    @assert occursin("n_hamiltonian_matrix_size", line)
    nnz = parse(Int64, split(line)[end]) - 1

    # Subtracting 1 to ignore the invalid
    # [999999999  999999999  999999999] cell
    line = readline(f)
    @assert occursin("n_cells_in_hamiltonian", line)
    n_cells = parse(Int64, split(line)[end]) - 1

    line = readline(f)
    @assert occursin("n_basis", line)
    n_basis = parse(Int64, split(line)[end])

    cells = SVector{3, Int}[]
    col_cell_ptr = Array{Int, 3}(undef, 2, n_cells, n_basis)
    row_ind = Vector{Int}(undef, nnz)

    @assert occursin("cell_index", readline(f))
    for _ in 1:n_cells
        push!(cells, SVector{3}(parse.(Int, split(readline(f)))))
    end

    # Jump over the invalid cell
    readline(f)

    @assert occursin("index_hamiltonian(1,:,:)", readline(f))
    for i in 1:n_cells
        col_cell_ptr[1, i, :] = parse.(Int, split(readline(f)))
    end

    # Jump over the invalid cell
    readline(f)

    @assert occursin("index_hamiltonian(2,:,:)", readline(f))
    for i in 1:n_cells
        col_cell_ptr[2, i, :] = parse.(Int, split(readline(f)))
    end

    # Jump over the invalid cell
    readline(f)

    @assert occursin("column_index_hamiltonian", readline(f))
    for i in 1:nnz
        row_ind[i] = parse(Int64, strip(readline(f)))
    end

    close(f)

    return FHIaimsOperatorMetadata(n_basis, row_ind, col_cell_ptr, cells)
end

function load_operator_data(dir::AbstractString, operatorkind::AbstractOperatorKind, ::Type{FHIaimsOperator})
    @debug "Loading operator data"

    ps = joinpath.(dir, filenames(operatorkind, Val(operatorkind.tag), FHIaimsOperator))
    ps = ps[ispath.(ps)]
    @argcheck !isempty(ps)

    exts = getindex.(splitext.(ps), 2)

    # TODO: might want to move those to a separate function as well
    # to allow loading an operator with (path, ...) instead of (dir, operatorkind, ...)
    for (p, ext) in zip(ps, exts)
        if ext == ".h5"
            io = h5open(p, "r")
            data = read(io["sparse_matrix"])
            close(io)

            # Ignore the last element (which is always zero)
            pop!(data)
            return data
        elseif ext == ".out"
            data = readdlm(p)

            # Reshape from (N, 1) to (N,)
            # and ignore the last element (which is always zero)
            data = reshape(data, :)
            pop!(data)
            return data
        else
            throw(error("Reading the file $(basename(p)) is unsupported"))
        end
    end
end

function FHIaimsOperator(dir::AbstractString, operatorkind::AbstractOperatorKind, metadata::FHIaimsOperatorMetadata)
    data = load_operator_data(dir, operatorkind, FHIaimsOperator)
    return FHIaimsOperator(operatorkind, data, metadata)
end

function load_operators(dir::AbstractString, operatorkinds, format::Type{FHIaimsOperator})
    metadata = load_operator_metadata(dir, format)
    operators = Dict(kind => FHIaimsOperator(dir, kind, metadata) for kind in operatorkinds)
    return operators, metadata
end

# Available operatorkinds of the input format
avail_operatorkinds(::Type{FHIaimsOperator}) = [Hamiltonian(:ref), Overlap(:ref)]

# Possible filenames for a given operatorkind and format
filenames(::Hamiltonian, ::Val{:ref}, ::Type{FHIaimsOperator}) = ["rs_hamiltonian.h5", "rs_hamiltonian.out"]
filenames(::Overlap, ::Val{:ref}, ::Type{FHIaimsOperator}) = ["rs_overlap.h5", "rs_overlap.out"]

# Find available operators in the supplied dir
function find_operatorkinds(dir::AbstractString, format::Type{<:AbstractOperator})
    found_operatorkinds = AbstractOperatorKind[]
    for operatorkind in avail_operatorkinds(format)
        names = filenames(operatorkind, Val(operatorkind.tag), format)
        ps = joinpath.(dir, names)
        if any(ispath.(ps))
            push!(found_operatorkinds, operatorkind)
        end
    end
    return found_operatorkinds
end

end