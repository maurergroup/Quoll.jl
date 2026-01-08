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

function build_operator_data(
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

### WRITING ###

function write_operators(
    ::Type{M}, dir::AbstractString, operators
) where {M<:DeepHBlockRealMetadata}
    ispath(dir) || mkpath(dir)

    # Metadata is shared -> write only once
    write_operator_metadata(dir, first(operators))

    for operator in operators
        write_operator_data(dir, operator)
    end
end

function write_operator_metadata(
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

function write_operator_data(
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
