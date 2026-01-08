struct FHIaimsSource <: AbstractSource end

### SPHERICAL HARMONICS CONVENTION ###

#! format: off
# Same order as wiki
# For m > 0 opposite CS phase
@memoize function default_shconv(source::FHIaimsSource)
    return SHConvention(
        [
                            [1],
                        [1,  2,  3],
                    [1,  2,  3,  4,  5], 
                [1,  2,  3,  4,  5,  6,  7],
            [1,  2,  3,  4,  5,  6,  7,  8,  9],
        ],
        [
                            [1], 
                        [1,  1, -1],
                    [1,  1,  1, -1,  1],
                [1,  1,  1,  1, -1,  1, -1],
            [1,  1,  1,  1,  1, -1,  1, -1,  1],
        ],
    )
end
#! format: on

### METADATA ALIASES ###

const FHIaimsCSCNoSpinRealMetadata{
O<:OperatorKind,
X<:FHIaimsSource,
S<:CSCRealSparsity,
B<:BasisSetMetadata,
Y<:SHConvention,
A<:AbstractSystem
} = RealMetadata{O,X,S,B,Y,A}

const FHIaimsCSCSpinRealMetadata{
O<:OperatorKind,
X<:FHIaimsSource,
S<:CSCRealSparsity,
B<:BasisSetMetadata,
Y<:SHConvention,
A<:AbstractSystem,
P<:SpinsMetadata
} = SpinRealMetadata{O,X,S,B,Y,A,P}

const FHIaimsCSCRealMetadata{
O<:OperatorKind,
X<:FHIaimsSource,
S<:CSCRealSparsity,
B<:BasisSetMetadata,
Y<:SHConvention,
A<:AbstractSystem,
P<:SpinsMetadata
} = Union{
    <:FHIaimsCSCNoSpinRealMetadata{O,X,S,B,Y,A},
    <:FHIaimsCSCSpinRealMetadata{O,X,S,B,Y,A,P},
}
op_data_type(::Type{<:FHIaimsCSCRealMetadata}) = CSCRealData
op_source_type(::Type{<:FHIaimsCSCRealMetadata}) = FHIaimsSource
op_sparsity_type(::Type{<:FHIaimsCSCRealMetadata}) = CSCRealSparsity

const FHIaimsDenseRecipMetadata{
O<:OperatorKind,
X<:FHIaimsSource,
S<:DenseRecipSparsity,
B<:BasisSetMetadata,
Y<:SHConvention,
A<:AbstractSystem,
K<:SVector{3}
} = RecipMetadata{O,X,S,B,Y,A,K}
op_data_type(::Type{<:FHIaimsDenseRecipMetadata}) = DenseRecipData
op_source_type(::Type{<:FHIaimsDenseRecipMetadata}) = FHIaimsSource
op_sparsity_type(::Type{<:FHIaimsDenseRecipMetadata}) = DenseRecipSparsity

### LOADING ###

function load_operator_basic_metadata(
    ::Type{<:FHIaimsCSCRealMetadata}, dir::AbstractString, kind::OperatorKind
)
    source = FHIaimsSource()
    atoms = load_atoms(source, dir)
    sparsity = CSCRealSparsity(source, dir)
    basisset = BasisSetMetadata(source, dir, atoms)
    shconv = default_shconv(source)
    return BasicMetadataContainer(kind, source, sparsity, basisset, shconv, atoms)
end

function load_operator_metadata(
    ::Type{<:FHIaimsCSCNoSpinRealMetadata},
    ::AbstractString,
    basic_metadata::BasicMetadataContainer,
)
    return RealMetadata(basic_metadata)
end

function load_operator_metadata(
    ::Type{<:FHIaimsCSCSpinRealMetadata},
    ::AbstractString,
    basic_metadata::BasicMetadataContainer,
)
    kind = op_kind(basic_metadata)
    source = op_source(basic_metadata)
    basisset = op_basisset(basic_metadata)

    spins = SpinsMetadata(source, kind, basisset)
    return SpinRealMetadata(basic_metadata, spins)
end

function load_operator_data(
    ::Type{M},
    dir::AbstractString,
    operatorkind::OperatorKind{K},
) where {M<:FHIaimsCSCRealMetadata,K}
    ps = joinpath.(dir, get_avail_filenames(M, operatorkind))
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

    return wrap_data(M, data)
end

@memoize function load_atoms(source::FHIaimsSource, dir::AbstractString)
    path = joinpath(dir, "geometry.in")
    atoms = load_system(path)

    # Recentre atom positions to be consistent with relative positions
    # that are used during real-space matrix integration in FHI-aims
    all(periodicity(atoms)) && (atoms = recentre(atoms))

    return atoms
end

@memoize function BasisSetMetadata(
    source::FHIaimsSource, dir::AbstractString, atoms::AbstractSystem
)
    p = joinpath(dir, "basis-indices.out")
    @argcheck ispath(p)

    atom2species = species(atoms, :)
    basis = Base.ImmutableDict(
        (
            species => BasisMetadata{Base.ImmutableDict{Symbol,Symbol}}[]
            for species in unique(atom2species)
        )...,
    )

    f = open(p, "r")

    # Reads both the empty lines and the header
    while isempty(strip(readline(f)))
    end

    species2firstatom = dictionary(
        species => findfirst(x -> x == species, atom2species)
        for species in unique(atom2species)
    )

    ib = 1
    while !(eof(f))
        _, type, iat, n, l, m = convert.(String, split(readline(f)))
        iat, n, l, m = parse.(Int, (iat, n, l, m))
        type = Symbol(type)

        if iat == species2firstatom[atom2species[iat]]
            push!(
                basis[atom2species[iat]],
                BasisMetadata(
                    atom2species[iat], n, l, m, Base.ImmutableDict(:type => type)
                ),
            )
        end

        ib += 1
    end
    close(f)

    # Lift arbitrary degeneracies if any
    basis = lift_arbitrary_degeneracy(basis)

    out_basisset = BasisSetMetadata(basis, atom2species)

    # Reorder basis functions according to FHIaims spherical harmonics convention
    out_basisset = convert_basisset_shconv(out_basisset, default_shconv(source))

    return out_basisset
end

@memoize function CSCRealSparsity(source::FHIaimsSource, dir::AbstractString)
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

    cells = SVector{3,Int}[]
    colcellptr = Array{Int,3}(undef, 2, n_cells, n_basis)
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

    return CSCRealSparsity(rowval, colcellptr, cells, true)
end

### PARSING ###

get_readformat(::Val{:fhiaims}) = FHIaimsCSCNoSpinRealMetadata
get_readformat(::Val{:fhiaimsspin}) = FHIaimsCSCSpinRealMetadata

get_avail_operatorkinds(::Type{<:FHIaimsCSCNoSpinRealMetadata}) = [
    Hamiltonian(; source=:ref),
    Overlap(; source=:ref),
]

get_avail_operatorkinds(::Type{<:FHIaimsCSCSpinRealMetadata}) = [
    Hamiltonian(; source=:ref, spin=:up),
    Hamiltonian(; source=:ref, spin=:down),
    Overlap(; source=:ref, spin=:up),
    Overlap(; source=:ref, spin=:down),
]

function get_avail_filenames(
    ::Type{M}, operatorkind::OperatorKind
) where {M<:FHIaimsCSCRealMetadata}
    return [
        get_avail_filename(M, statictuple(operatorkind), Val(ext))
        for ext in (:hdf5, :plain)
    ]
end

# Can attach additional traits to the unionall constant and dispatch using holy traits
# FHIaimsCSCNoSpinRealOperator

get_avail_filename(
    ::Type{<:FHIaimsCSCNoSpinRealMetadata},
    ::Tuple{Val{:Hamiltonian},Pair{Val{:source},Val{:ref}}},
    ::Val{:hdf5},
) = "rs_hamiltonian.h5"

get_avail_filename(
    ::Type{<:FHIaimsCSCNoSpinRealMetadata},
    ::Tuple{Val{:Hamiltonian},Pair{Val{:source},Val{:ref}}},
    ::Val{:plain},
) = "rs_hamiltonian.out"

get_avail_filename(
    ::Type{<:FHIaimsCSCNoSpinRealMetadata},
    ::Tuple{Val{:Overlap},Pair{Val{:source},Val{:ref}}},
    ::Val{:hdf5},
) = "rs_overlap.h5"

get_avail_filename(
    ::Type{<:FHIaimsCSCNoSpinRealMetadata},
    ::Tuple{Val{:Overlap},Pair{Val{:source},Val{:ref}}},
    ::Val{:plain},
) = "rs_overlap.out"

# FHIaimsCSCSpinRealOperator

get_avail_filename(
    ::Type{<:FHIaimsCSCSpinRealMetadata},
    ::Tuple{Val{:Hamiltonian},Pair{Val{:source},Val{:ref}},Pair{Val{:spin},Val{:up}}},
    ::Val{:hdf5},
) = "rs_hamiltonian_up.h5"

get_avail_filename(
    ::Type{<:FHIaimsCSCSpinRealMetadata},
    ::Tuple{Val{:Hamiltonian},Pair{Val{:source},Val{:ref}},Pair{Val{:spin},Val{:up}}},
    ::Val{:plain},
) = "rs_hamiltonian_up.out"

get_avail_filename(
    ::Type{<:FHIaimsCSCSpinRealMetadata},
    ::Tuple{Val{:Hamiltonian},Pair{Val{:source},Val{:ref}},Pair{Val{:spin},Val{:down}}},
    ::Val{:hdf5},
) = "rs_hamiltonian_down.h5"

get_avail_filename(
    ::Type{<:FHIaimsCSCSpinRealMetadata},
    ::Tuple{Val{:Hamiltonian},Pair{Val{:source},Val{:ref}},Pair{Val{:spin},Val{:down}}},
    ::Val{:plain},
) = "rs_hamiltonian_dn.out"

get_avail_filename(
    ::Type{<:FHIaimsCSCSpinRealMetadata},
    ::Tuple{Val{:Overlap},Pair{Val{:source},Val{:ref}},Pair{Val{:spin},S}},
    ::Val{:hdf5},
) where {S<:Union{Val{:up},Val{:down}}} = "rs_overlap.h5"

get_avail_filename(
    ::Type{<:FHIaimsCSCSpinRealMetadata},
    ::Tuple{Val{:Overlap},Pair{Val{:source},Val{:ref}},Pair{Val{:spin},S}},
    ::Val{:plain},
) where {S<:Union{Val{:up},Val{:down}}} = "rs_overlap.out"
