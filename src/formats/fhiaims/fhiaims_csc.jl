using DelimitedFiles: readdlm

struct FHIaimsCSCMetadata{E} <: AbstractFHIaimsMetadata
    sparsity::RealCSCSparsity
    basisset::BasisSetMetadata{E}
end

struct FHIaimsCSCOperator{O<:AbstractOperatorKind, T<:AbstractFloat, E} <: AbstractFHIaimsOperator
    kind::O
    data::Vector{T}
    metadata::FHIaimsCSCMetadata{E}
end

get_readformat(::Val{:fhiaims}) = FHIaimsCSCOperator

get_avail_operatorkinds(::Type{FHIaimsCSCOperator}) = [Hamiltonian(:ref), Overlap(:ref)]
get_avail_filenames(::Hamiltonian, ::Val{:ref}, ::Type{FHIaimsCSCOperator}) = ["rs_hamiltonian.h5", "rs_hamiltonian.out"]
get_avail_filenames(::Overlap, ::Val{:ref}, ::Type{FHIaimsCSCOperator}) = ["rs_overlap.h5", "rs_overlap.out"]

function RealCSCSparsity(dir::AbstractString, ::Type{FHIaimsCSCOperator})
    @debug "Loading operator sparsity"

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
    colcellptr = Array{Int, 3}(undef, 2, n_cells, n_basis)
    rowval = Vector{Int}(undef, nnz)

    @assert occursin("cell_index", readline(f))
    for _ in 1:n_cells
        push!(cells, SVector{3}(parse.(Int, split(readline(f)))))
    end

    # Jump over the invalid cell
    readline(f)

    @assert occursin("index_hamiltonian(1,:,:)", readline(f))
    for i in 1:n_cells
        colcellptr[1, i, :] = parse.(Int, split(readline(f)))
    end

    # Jump over the invalid cell
    readline(f)

    @assert occursin("index_hamiltonian(2,:,:)", readline(f))
    for i in 1:n_cells
        colcellptr[2, i, :] = parse.(Int, split(readline(f)))
    end

    # Jump over the invalid cell
    readline(f)

    @assert occursin("column_index_hamiltonian", readline(f))
    for i in 1:nnz
        rowval[i] = parse(Int64, strip(readline(f)))
    end

    close(f)

    return RealCSCSparsity(rowval, colcellptr, cells)
end

function load_operators(dir::AbstractString, atoms::AbstractSystem, operatorkinds, ::Type{FHIaimsCSCOperator})
    metadata = FHIaimsCSCMetadata(dir, atoms)
    operators = Dict(kind => FHIaimsCSCOperator(dir, kind, metadata) for kind in operatorkinds)
    return operators
end

# TODO: Not sure if I'll keep this function, it could be useful to have if
# we want to use load_operator_metadata in quoll/bin
function load_operator_metadata(dir::AbstractString, atoms::AbstractSystem, ::Type{FHIaimsCSCOperator})
    return FHIaimsCSCMetadata(dir, atoms)
end

function FHIaimsCSCMetadata(dir::AbstractString, atoms::AbstractSystem)
    @debug "Loading operator metadata"
    sparsity = RealCSCSparsity(dir, FHIaimsCSCOperator)
    basisset = BasisSetMetadata(dir, atoms, FHIaimsCSCOperator)
    return FHIaimsCSCMetadata(sparsity, basisset)
end

function FHIaimsCSCOperator(dir::AbstractString, operatorkind::AbstractOperatorKind, metadata::FHIaimsCSCMetadata)
    data = load_operator_data(dir, operatorkind, FHIaimsCSCOperator)
    return FHIaimsCSCOperator(operatorkind, data, metadata)
end

function load_operator_data(dir::AbstractString, operatorkind::AbstractOperatorKind, ::Type{FHIaimsCSCOperator})
    @debug "Loading operator data"

    ps = joinpath.(dir, get_avail_filenames(operatorkind, Val(operatorkind.tag), FHIaimsCSCOperator))
    ps = ps[ispath.(ps)]
    @argcheck !isempty(ps)

    exts = getindex.(splitext.(ps), 2)

    # Choose the first filename
    p = ps[1]
    ext = exts[1]

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
