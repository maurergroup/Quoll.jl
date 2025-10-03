# TODO: should probably define an API to make it easier for developers
# to implement new formats 
module OperatorIO
using StaticArrays

export read_format, write_format, operator_kind
export AbstractOperator
export FHIaimsOperator, DeepHOperator
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
struct RealBSparseOperator <: AbstractBSparseOperator end

struct FHIaimsOperatorMetadata <: AbstractOperatorMetadata
    row_ind::Vector{Int}
    col_cell_ptr::Array{Int, 3}
    cells::Vector{SVector{3, Int}}
end

# This will not be a singleton, which means if we want to dispatch
# on the format without passing values into a method we will need to dispatch
# on the type instead of an instance (::Type{FHIaimsOperator}).
# Would that be inconsistent if I use singleton instances for e.g. FC99V?
struct FHIaimsOperator{O<:AbstractOperatorKind, T<:AbstractFloat} <: AbstractOperator
    operator::O
    data::Vector{T}
    metadata::FHIaimsOperatorMetadata
end

struct DeepHOperator <: AbstractOperator end

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
    
    open(joinpath(dir, "rs_indices.out"), "r") do f
        line = readline(f)
        @assert occursin("n_hamiltonian_matrix_size", readline(f))
        # Subtracting 1 to ignore one element (last element)
        # which is always zero
        nnz = parse(Int64, split(line)[end]) - 1

        line = readline(f)
        @assert occursin("n_cells_in_hamiltonian", line)
        # Subtracting 1 to ignore the invalid
        # [999999999  999999999  999999999] cell
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

    end
    
end

function FHIaimsOperator()
end

function load_operators()
end

end