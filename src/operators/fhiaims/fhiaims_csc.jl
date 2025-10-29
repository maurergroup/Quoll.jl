using DelimitedFiles: readdlm

# TODO: What about atoms, should they be attached to the operator metadata?
# Some atom information is already present in metadata,
# e.g. BasisSetMetadata contains chemical species and atom indices.
# For completeness it might make sense to attach atoms fully.
# On the other hand, what if someone wants to initialise the operator
# without information on atoms? Is that even a concern?
# Also, type of atoms would make the operator type even more parametric,
# which could lead to longer compilation times if typeof(atoms) changes.
# However, that's probably unlikely because although typeof(atoms) is
# itself a parametric type of subtype AbstractSystem, the parameters define
# the dimension of space and Unitful quantities which, for example, would not
# change for different atoms in a dataset

struct FHIaimsCSCMetadata{A<:AbstractSystem, E} <: AbstractFHIaimsMetadata
    atoms::A
    sparsity::RealCSCSparsity
    basisset::BasisSetMetadata{E}
    spins::Union{SpinsMetadata, Nothing}
end

struct FHIaimsCSCOperator{O<:OperatorKind, T<:AbstractFloat, A<:AbstractSystem, E} <: AbstractFHIaimsOperator
    kind::O
    data::Vector{T}
    metadata::FHIaimsCSCMetadata{A, E}
end

get_float(operator::FHIaimsCSCOperator) = typeof(operator).parameters[2]

get_readformat(::Val{:fhiaims}) = FHIaimsCSCOperator

get_avail_operatorkinds(::Type{FHIaimsCSCOperator}) = [
    Hamiltonian(source = :ref, spin = :none),
    Hamiltonian(source = :ref, spin = :up),
    Hamiltonian(source = :ref, spin = :down),
    Overlap(source = :ref),
]

get_avail_filenames(
    ::Hamiltonian,
    ::Pair{Val{:source}, Val{:ref}},
    ::Pair{Val{:spin}, Val{:none}},
    ::Type{FHIaimsCSCOperator}
) = ["rs_hamiltonian.h5", "rs_hamiltonian.out"]

get_avail_filenames(
    ::Hamiltonian,
    ::Pair{Val{:source}, Val{:ref}},
    ::Pair{Val{:spin}, Val{:up}},
    ::Type{FHIaimsCSCOperator}
) = ["rs_hamiltonian_up.h5", "rs_hamiltonian_up.out"]

get_avail_filenames(
    ::Hamiltonian,
    ::Pair{Val{:source}, Val{:ref}},
    ::Pair{Val{:spin}, Val{:down}},
    ::Type{FHIaimsCSCOperator}
) = ["rs_hamiltonian_down.h5", "rs_hamiltonian_dn.out"]

get_avail_filenames(
    ::Overlap,
    ::Pair{Val{:source}, Val{:ref}},
    ::Type{FHIaimsCSCOperator}
) = ["rs_overlap.h5", "rs_overlap.out"]

function RealCSCSparsity(dir::AbstractString, ::Type{FHIaimsCSCOperator})
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

    return RealCSCSparsity(rowval, colcellptr, cells, true)
end

function convert_sparsity(in_metadata::FHIaimsCSCMetadata, out_sparsity_type::Type{RealBlockSparsity}; hermitian = true)
    return convert_sparsity(get_sparsity(in_metadata), get_basisset(in_metadata), out_sparsity_type, hermitian = hermitian)
end

function FHIaimsCSCOperator(dir::AbstractString, kind::OperatorKind)
    return FHIaimsCSCOperator(
        kind,
        load_operator_data(dir, kind, FHIaimsCSCOperator),
        load_operator_metadata(dir, kind, FHIaimsCSCOperator)
    )
end

function load_operator_data(dir::AbstractString, operatorkind::OperatorKind, ::Type{FHIaimsCSCOperator})
    ps = joinpath.(dir, get_avail_filenames(operatorkind, FHIaimsCSCOperator))
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

function load_operator_metadata(dir::AbstractString, kind::OperatorKind, ::Type{FHIaimsCSCOperator})
    return FHIaimsCSCMetadata(dir, kind) 
end

function FHIaimsCSCMetadata(dir::AbstractString, kind::OperatorKind)
    atoms = load_atoms(dir, FHIaimsCSCOperator)
    sparsity = RealCSCSparsity(dir, FHIaimsCSCOperator)
    basisset = BasisSetMetadata(dir, atoms, FHIaimsCSCOperator)
    spins = SpinsMetadata(kind, basisset, FHIaimsCSCOperator)
    return FHIaimsCSCMetadata(atoms, sparsity, basisset, spins)
end

# Specialised load_operators for loading metadata once because we know that
# all FHIaimsCSCOperators in a dir will share metadata to some extent
function load_operators(dir::AbstractString, operatorkinds, ::Type{FHIaimsCSCOperator})
    # Load shared metadata
    atoms = load_atoms(dir, FHIaimsCSCOperator)
    sparsity = RealCSCSparsity(dir, FHIaimsCSCOperator)
    basisset = BasisSetMetadata(dir, atoms, FHIaimsCSCOperator)

    operators = []
    for kind in operatorkinds
        data = load_operator_data(dir, kind, FHIaimsCSCOperator)
        spinsset = SpinsMetadata(kind, basisset, FHIaimsCSCOperator)
        metadata = FHIaimsCSCMetadata(atoms, sparsity, basisset, spinsset)
        push!(operators, FHIaimsCSCOperator(kind, data, metadata))
    end

    return operators
end
