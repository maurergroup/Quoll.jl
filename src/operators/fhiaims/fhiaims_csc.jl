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

struct FHIaimsCSCMetadata{
    A<:AbstractSystem,
    S<:RealCSCSparsity,
    B<:BasisSetMetadata,
    P<:Union{SpinsMetadata, Nothing}
} <: AbstractFHIaimsMetadata{A, S, B, P}
    atoms::A
    sparsity::S
    basisset::B
    spins::P
end

struct FHIaimsCSCOperator{
    O<:OperatorKind,
    T<:Number,
    D<:AbstractVector{T},
    M<:FHIaimsCSCMetadata
} <: AbstractFHIaimsOperator{O, T, D, M}
    kind::O
    data::D
    metadata::M
end

### READING METHODS ###

function FHIaimsCSCOperator(dir::AbstractString, kind::OperatorKind)
    return FHIaimsCSCOperator(
        kind,
        load_operator_data(dir, kind, FHIaimsCSCOperator),
        load_operator_metadata(dir, kind, FHIaimsCSCOperator)
    )
end

function load_operator_data(dir::AbstractString, operatorkind::OperatorKind{K}, ::Type{FHIaimsCSCOperator}) where K
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
    elseif ext == ".out"
        data = readdlm(p)

        # Reshape from (N, 1) to (N,)
        # and ignore the last element (which is always zero)
        data = reshape(data, :)
        pop!(data)
    else
        throw(error("Reading the file $(basename(p)) is unsupported"))
    end

    if K == :Hamiltonian
        data *= HARTREE_CODATA_2002
    end

    return data
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

### APP-RELATED METHODS ###

get_readformat(::Val{:fhiaims}) = FHIaimsCSCOperator

get_avail_operatorkinds(::Type{FHIaimsCSCOperator}) = [
    Hamiltonian(source = :ref, spin = :none),
    Hamiltonian(source = :ref, spin = :up),
    Hamiltonian(source = :ref, spin = :down),
    Overlap(source = :ref),
]

function get_avail_filenames(operatorkind::OperatorKind, pairtypes::Tuple, T::Type{FHIaimsCSCOperator})
    return [get_avail_filename(operatorkind, pairtypes..., Val(output), T) for output in (:hdf5, :plain)]
end

get_avail_filename(
    ::Hamiltonian,
    ::Pair{Val{:source}, Val{:ref}},
    ::Pair{Val{:spin}, Val{:none}},
    ::Val{:hdf5},
    ::Type{FHIaimsCSCOperator}
) = "rs_hamiltonian.h5"

get_avail_filename(
    ::Hamiltonian,
    ::Pair{Val{:source}, Val{:ref}},
    ::Pair{Val{:spin}, Val{:none}},
    ::Val{:plain},
    ::Type{FHIaimsCSCOperator}
) = "rs_hamiltonian.out"

get_avail_filename(
    ::Hamiltonian,
    ::Pair{Val{:source}, Val{:ref}},
    ::Pair{Val{:spin}, Val{:up}},
    ::Val{:hdf5},
    ::Type{FHIaimsCSCOperator}
) = "rs_hamiltonian_up.h5"

get_avail_filename(
    ::Hamiltonian,
    ::Pair{Val{:source}, Val{:ref}},
    ::Pair{Val{:spin}, Val{:up}},
    ::Val{:plain},
    ::Type{FHIaimsCSCOperator}
) = "rs_hamiltonian_up.out"

get_avail_filename(
    ::Hamiltonian,
    ::Pair{Val{:source}, Val{:ref}},
    ::Pair{Val{:spin}, Val{:down}},
    ::Val{:hdf5},
    ::Type{FHIaimsCSCOperator}
) = "rs_hamiltonian_down.h5"

get_avail_filename(
    ::Hamiltonian,
    ::Pair{Val{:source}, Val{:ref}},
    ::Pair{Val{:spin}, Val{:down}},
    ::Val{:plain},
    ::Type{FHIaimsCSCOperator}
) = "rs_hamiltonian_dn.out"

get_avail_filename(
    ::Overlap,
    ::Pair{Val{:source}, Val{:ref}},
    ::Val{:hdf5},
    ::Type{FHIaimsCSCOperator}
) = "rs_overlap.h5"

get_avail_filename(
    ::Overlap,
    ::Pair{Val{:source}, Val{:ref}},
    ::Val{:plain},
    ::Type{FHIaimsCSCOperator}
) = "rs_overlap.out"
