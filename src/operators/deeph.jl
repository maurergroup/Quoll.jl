struct DeepHSource <: AbstractSource end

#! format: off
@memoize function default_shconv(source::DeepHSource)
    return SHConvention(
        [
                     [1],
                  [3, 1, 2],
               [3, 5, 1, 4, 2], 
            [4, 5, 3, 6, 2, 7, 1],
        ],
        [
                     [1],
                  [1, 1, 1],
               [1, 1, 1, 1, 1],
            [1, 1, 1, 1, 1, 1, 1],
        ],
    )
end
#! format: on

### DATA ALIASES ###

const DeepHBlockRealData{T} = DataContainer{
    T,2,
    <:AbstractDictionary{NTuple{5,Int},<:AbstractMatrix{T}},
    <:DeepHSource,
    <:BlockRealSparsity,
}

### METADATA ALIASES ###

# DeepHBlockRealMetadata

const DeepHBlockNoSpinRealMetadata{
O<:OperatorKind,
X<:DeepHSource,
S<:BlockRealSparsity,
B<:BasisSetMetadata,
Y<:SHConvention,
A<:AbstractSystem
} = RealMetadata{O,X,S,B,Y,A}

const DeepHBlockSpinRealMetadata{
O<:OperatorKind,
X<:DeepHSource,
S<:BlockRealSparsity,
B<:BasisSetMetadata,
Y<:SHConvention,
A<:AbstractSystem,
P<:SpinsMetadata
} = SpinRealMetadata{O,X,S,B,Y,A,P}

const DeepHBlockRealMetadata{
O<:OperatorKind,
X<:DeepHSource,
S<:BlockRealSparsity,
B<:BasisSetMetadata,
Y<:SHConvention,
A<:AbstractSystem,
P<:SpinsMetadata
} = Union{
    <:DeepHBlockNoSpinRealMetadata{O,X,S,B,Y,A},
    <:DeepHBlockSpinRealMetadata{O,X,S,B,Y,A,P},
}
op_data_type(::Type{<:DeepHBlockRealMetadata}) = DeepHBlockRealData
op_source_type(::Type{<:DeepHBlockRealMetadata}) = DeepHSource
op_sparsity_type(::Type{<:DeepHBlockRealMetadata}) = BlockRealSparsity

### FACTORIES ###

function build_data(
    ::Type{<:DeepHBlockRealData}, metadata::M, value::T, initialised::Bool
) where {M<:AbstractMetadata,T<:Number}
    sparsity = op_sparsity(metadata)
    deeph_keys = get_deeph_keys(sparsity)

    basisset = op_basisset(metadata)
    atom2nbasis = get_atom2nbasis(basisset)

    data = map(Indices(deeph_keys)) do k
        _, _, _, iat, jat = k
        array = Matrix{T}(undef, atom2nbasis[iat], atom2nbasis[jat])
        initialised && fill!(array, value)
        return array
    end

    return wrap_data(M, data)
end

### LOADING ###

function load_metadata_basic(
    ::Type{M}, dir::AbstractString, kind::OperatorKind
) where {M<:DeepHBlockRealMetadata}
    source = DeepHSource()
    atoms = load_atoms(source, dir)
    sparsity = BlockRealSparsity(M, dir, kind)
    basisset = BasisSetMetadata(source, dir, atoms)
    shconv = default_shconv(source)
    return BasicMetadataContainer(kind, source, sparsity, basisset, shconv, atoms)
end

function load_metadata(
    ::Type{<:DeepHBlockNoSpinRealMetadata},
    ::AbstractString,
    basic_metadata::BasicMetadataContainer,
)
    return RealMetadata(basic_metadata)
end

function load_metadata(
    ::Type{<:DeepHBlockSpinRealMetadata},
    ::AbstractString,
    basic_metadata::BasicMetadataContainer,
)
    kind = op_kind(basic_metadata)
    source = op_source(basic_metadata)
    basisset = op_basisset(basic_metadata)

    spins = SpinsMetadata(source, kind, basisset)
    return SpinRealMetadata(basic_metadata, spins)
end

function load_data(
    ::Type{M},
    dir::AbstractString,
    operatorkind::OperatorKind,
) where {M<:DeepHBlockRealMetadata}

    p = joinpath(dir, get_avail_filename(M, statictuple(operatorkind)))
    @argcheck ispath(p)

    io = h5open(p, "r")
    V = typeof(read(first(io)))

    return load_data(M, io, V)
end

function load_data(
    ::Type{M}, io::HDF5.File, ::Type{V}
) where {M<:DeepHBlockRealMetadata,V<:AbstractMatrix{T}} where {T}
    keys_list = NTuple{5,Int}[]
    values_list = V[]

    for (key_str, value) in pairs(io)
        key = JSON.parse(key_str, NTuple{5,Int})
        block = read(value, T)
        push!(keys_list, key)
        push!(values_list, permutedims(block, (2, 1)))
    end
    close(io)

    data = Dictionary(keys_list, values_list)
    
    return wrap_data(M, data)
end

@memoize function load_atoms(source::DeepHSource, dir::AbstractString)
    elements = readdlm(joinpath(dir, "element.dat"), Int)

    # 3xN
    positions = readdlm(joinpath(dir, "site_positions.dat"))

    # 3x3 with columns as lattice vectors
    lattice = readdlm(joinpath(dir, "lat.dat"))

    atoms = periodic_system(
        [
            element => position * u"Å"
            for (element, position) in zip(elements, eachcol(positions))
        ],
        [cell_vector * u"Å" for cell_vector in eachcol(lattice)],
    )

    return atoms
end

@memoize function BasisSetMetadata(
    source::DeepHSource, dir::AbstractString, atoms::AbstractSystem
)
    p = joinpath(dir, "orbital_types.dat")
    @argcheck ispath(p)

    atom2species = species(atoms, :)
    basis = Base.ImmutableDict(
        (
            species => BasisMetadata{Nothing}[]
            for species in unique(atom2species)
        )...,
    )

    species2firstatom = dictionary(
        species => findfirst(x -> x == species, atom2species)
        for species in unique(atom2species)
    )

    for (iat, line) in enumerate(readlines(p))
        z = atom2species[iat]
        if iat == species2firstatom[z]
            l_counts = Dict{Int,Int}()
            for l in parse.(Int, split(line))
                l ∈ keys(l_counts) || push!(l_counts, l => 1)
                n = l_counts[l]
                for m in ((-l):l)
                    push!(
                        basis[z],
                        BasisMetadata(
                            z, n, l, m, nothing
                        ),
                    )
                end
                l_counts[l] += 1
            end
        end
    end

    out_basisset = BasisSetMetadata(basis, atom2species)

    # Reorder basis functions according to DeepH spherical harmonics convention
    out_basisset = convert_basisset_shconv(out_basisset, default_shconv(source))

    return out_basisset
end

@memoize function BlockRealSparsity(
    format::Type{M}, dir::AbstractString, operatorkind::OperatorKind
) where {M<:DeepHBlockRealMetadata}
    p = joinpath(dir, get_avail_filename(M, statictuple(operatorkind)))
    @argcheck ispath(p)

    # Vector of "[Rx, Ry, Rz, i, j]" strings
    io = h5open(p, "r")
    keys_str = keys(io)
    close(io)

    images = SVector{3,Int}[]
    ij2images = Dictionary{NTuple{2,Int},Vector{SVector{3,Int}}}()
    for key_str in keys_str
        Rx, Ry, Rz, ij... = JSON.parse(key_str, NTuple{5,Int})
        R = SVector(Rx, Ry, Rz)

        get!(Vector{SVector{3,Int}}, ij2images, ij)
        push!(ij2images[ij], R)
        push!(images, R)
    end

    return BlockRealSparsity(ij2images, unique(images), false)
end

# load_metadata

# load_data

### WRITING ###

function write_operators(
    ::Type{M}, dir::AbstractString, operators
) where {M<:DeepHBlockRealMetadata}
    ispath(dir) || mkpath(dir)

    # Metadata is shared -> write only once
    write_metadata(dir, first(operators))

    for operator in operators
        write_data(dir, operator)
    end
end

function write_metadata(
    dir::AbstractString,
    operator::AbstractOperator{<:OperatorKind,<:DeepHBlockRealMetadata},
)
    write_rlist(dir, operator)
    write_element(dir, operator)
    write_info(dir, operator)
    write_lat(dir, operator)
    write_rlat(dir, operator)
    write_orbital_types(dir, operator)
    return write_site_positions(dir, operator)
end

function write_data(
    dir::AbstractString,
    operator::AbstractOperator{<:OperatorKind,<:DeepHBlockRealMetadata},
)
    metadata = op_metadata(operator)
    data_body = unwrap_data(op_data(operator))

    filename = first(get_avail_filenames(metadata))
    h5open(joinpath(dir, filename), "w") do db
        for (key, block) in pairs(data_body)
            key_str = string(collect(key))

            # Invert the layout of the array (from column major to row major)
            # This was done originally because this format is primarily used in Python
            write(db, key_str, permutedims(block, ndims(block):-1:1))
        end
    end
end

function write_rlist(
    dir::AbstractString,
    operator::AbstractOperator{<:OperatorKind,<:DeepHBlockRealMetadata},
)
    return writedlm(joinpath(dir, "R_list.dat"), op_images(op_sparsity(operator)))
end

function write_element(
    dir::AbstractString,
    operator::AbstractOperator{<:OperatorKind,<:DeepHBlockRealMetadata},
)
    atomic_numbers = atomic_number(op_atoms(operator), :)
    return writedlm(joinpath(dir, "element.dat"), atomic_numbers)
end

# TODO: will need an additional input which stores eigenvalues/eigenvectors
# (for fermi level)
function write_info(
    dir::AbstractString,
    operator::AbstractOperator{<:OperatorKind,<:DeepHBlockRealMetadata},
)
    kind = op_kind(operator)
    isspinful = haskey(kind.tags, :spin) && isequal(kind.spin, :soc)
    return JSON.json(
        joinpath(dir, "info.json"),
        Dict("isspinful" => isspinful);
        pretty=true,
    )
end

function write_lat(
    dir::AbstractString,
    operator::AbstractOperator{<:OperatorKind,<:DeepHBlockRealMetadata},
)
    lat = ustrip.(hcat(cell_vectors(op_atoms(operator))...))
    return writedlm(joinpath(dir, "lat.dat"), lat)
end

function write_rlat(
    dir::AbstractString,
    operator::AbstractOperator{<:OperatorKind,<:DeepHBlockRealMetadata},
)
    lat = ustrip.(hcat(cell_vectors(op_atoms(operator))...))
    rlat = transpose(lat) \ (2π * I(3))
    return writedlm(joinpath(dir, "rlat.dat"), rlat)
end

function write_orbital_types(
    dir::AbstractString,
    operator::AbstractOperator{<:OperatorKind,<:DeepHBlockRealMetadata},
)
    basisset = op_basisset(operator)
    open(joinpath(dir, "orbital_types.dat"), "w") do io
        for atom in op_atoms(operator)
            writedlm(
                io, transpose(get_angular_momenta(basis_species(basisset, species(atom))))
            )
        end
    end
end

function write_site_positions(
    dir::AbstractString,
    operator::AbstractOperator{<:OperatorKind,<:DeepHBlockRealMetadata},
)
    positions = ustrip.(hcat(position(op_atoms(operator), :)...))
    return writedlm(joinpath(dir, "site_positions.dat"), positions)
end

### MISC METHODS ###

function get_deeph_keys(sparsity::BlockRealSparsity)
    return [
        (Rx, Ry, Rz, iat, jat)
        for (iat, jat) in keys(sparsity.ij2images)
        for (Rx, Ry, Rz) in sparsity.ij2images[(iat, jat)]
    ]
end

### PARSING ###

get_readformat(::Val{:deeph}) = DeepHBlockNoSpinRealMetadata
get_readformat(::Val{:deephspin}) = DeepHBlockSpinRealMetadata

get_writeformat(::Val{:deeph}) = DeepHBlockNoSpinRealMetadata
get_writeformat(::Val{:deephspin}) = DeepHBlockSpinRealMetadata

get_avail_operatorkinds(::Type{<:DeepHBlockNoSpinRealMetadata}) = [
    Hamiltonian(; source=:ref),
    Hamiltonian(; source=:pred),
    Overlap(; source=:ref),
]

get_avail_operatorkinds(::Type{<:DeepHBlockSpinRealMetadata}) = [
    Hamiltonian(; source=:ref, spin=:soc),
    Hamiltonian(; source=:pred, spin=:soc),
    Overlap(; source=:ref, spin=:soc),
]

# DeepHBlockNoSpinRealMetadata

get_avail_filename(
    ::Type{<:DeepHBlockNoSpinRealMetadata},
    ::Tuple{Val{:Hamiltonian},Pair{Val{:source},Val{:ref}}},
) = "hamiltonians.h5"

get_avail_filename(
    ::Type{<:DeepHBlockNoSpinRealMetadata},
    ::Tuple{Val{:Hamiltonian},Pair{Val{:source},Val{:pred}}},
) = "hamiltonians_pred.h5"

get_avail_filename(
    ::Type{<:DeepHBlockNoSpinRealMetadata},
    ::Tuple{Val{:Overlap},Pair{Val{:source},Val{:ref}}},
) = "overlaps.h5"

# DeepHBlockSpinRealMetadata

get_avail_filename(
    ::Type{<:DeepHBlockSpinRealMetadata},
    ::Tuple{Val{:Hamiltonian},Pair{Val{:source},Val{:ref}},Pair{Val{:spin},Val{:soc}}},
) = "hamiltonians.h5"

get_avail_filename(
    ::Type{<:DeepHBlockSpinRealMetadata},
    ::Tuple{Val{:Hamiltonian},Pair{Val{:source},Val{:pred}},Pair{Val{:spin},Val{:soc}}},
) = "hamiltonians_pred.h5"

get_avail_filename(
    ::Type{<:DeepHBlockSpinRealMetadata},
    ::Tuple{Val{:Overlap},Pair{Val{:source},Val{:ref}},Pair{Val{:spin},Val{:soc}}},
) = "overlaps.h5"
