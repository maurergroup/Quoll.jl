const DeepHDict{T, N} = AbstractAnyDict{NTuple{5, Int}, <:AbstractArray{T, N}}

struct DeepHMetadata{A<:AbstractSystem, S<:RealBlockSparsity, B<:BasisSetMetadata, P<:Union{SpinsMetadata, Nothing}} <: AbstractOperatorMetadata{A, S, B, P}
    atoms::A
    sparsity::S
    basisset::B
    spins::P
end

struct DeepHOperator{O<:OperatorKind, T<:Number, D<:DeepHDict{T}, M<:DeepHMetadata} <: AbstractOperator{O, T, D, M}
    kind::O
    data::D
    metadata::M
end

const DeepHSHConversion = SHConversion(
    [[1], [3, 1, 2], [3, 5, 1, 4, 2], [4, 5, 3, 6, 2, 7, 1]],
    [[1], [1, 1, 1], [1, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1, 1]],
)

SHConversion(::Type{DeepHOperator}) = DeepHSHConversion

### WRITE METHODS ###

function write_operators(dir::AbstractString, operators::AbstractVector{<:DeepHOperator})
    ispath(dir) || mkpath(dir)

    # Metadata is shared -> write only once
    write_operator_metadata(dir, first(operators))

    for operator in operators
        write_operator_data(dir, operator)
    end
end

function write_operator_metadata(dir::AbstractString, operator::DeepHOperator)
    write_rlist(dir, operator)
    write_element(dir, operator)
    write_info(dir, operator)
    write_lat(dir, operator)
    write_rlat(dir, operator)
    write_orbital_types(dir, operator)
    write_site_positions(dir, operator)
end

function write_operator_data(dir, operator::DeepHOperator)
    h5open(joinpath(dir, get_avail_filename(operator)), "w") do db
        for (key, block) in pairs(get_data(operator))
            key_str = string(collect(key))

            # Invert the layout of the array (from column major to row major)
            write(db, key_str, permutedims(block, ndims(block):-1:1))
        end
    end
end

function write_rlist(dir::AbstractString, operator::DeepHOperator)
    writedlm(joinpath(dir, "R_list.dat"), get_sparsity(operator).images)
end

function write_element(dir::AbstractString, operator::DeepHOperator)
    atomic_numbers = atomic_number(get_atoms(operator), :)
    writedlm(joinpath(dir, "element.dat"), atomic_numbers)
end

# TODO: will need an additional input which stores eigenvalues/eigenvectors
# (for fermi level)
function write_info(dir::AbstractString, operator::DeepHOperator)
    JSON.json(
        joinpath(dir, "info.json"),
        Dict("isspinful" => get_kind(operator).spin == :soc);
        pretty = true
    )
end

function write_lat(dir::AbstractString, operator::DeepHOperator)
    lat = ustrip.(hcat(cell_vectors(get_atoms(operator))...))
    writedlm(joinpath(dir, "lat.dat"), lat)
end

# TODO: rlat will be needed elsewhere as well so move
# its construction into a separate function
function write_rlat(dir::AbstractString, operator::DeepHOperator)
    lat = ustrip.(hcat(cell_vectors(get_atoms(operator))...))
    rlat = transpose(lat) \ (2π * I(3))
    writedlm(joinpath(dir, "rlat.dat"), rlat)
end

function write_orbital_types(dir::AbstractString, operator::DeepHOperator)
    basisset = get_basisset(operator)
    open(joinpath(dir, "orbital_types.dat"), "w") do io
        for atom in get_atoms(operator)
            writedlm(io, transpose(get_angular_momenta(basis_species(basisset, species(atom)))))
        end
    end
end

function write_site_positions(dir::AbstractString, operator::DeepHOperator)
    positions = ustrip.(hcat(position(get_atoms(operator), :)...))
    writedlm(joinpath(dir, "site_positions.dat"), positions)
end

### APP-RELATED METHODS ###

get_writeformat(::Val{:deeph}) = DeepHOperator

get_avail_operatorkinds(::Type{DeepHOperator}) = [
    Hamiltonian(source = :ref, spin = :none),
    Hamiltonian(source = :ref, spin = :soc),
    Hamiltonian(source = :pred, spin = :none),
    Hamiltonian(source = :pred, spin = :soc),
    Overlap(source = :ref),
]

get_avail_filename(
    ::Hamiltonian,
    ::Pair{Val{:source}, Val{:ref}},
    ::Pair{Val{:spin}, Val{SPIN}},
    ::Type{<:DeepHOperator}
) where SPIN = "hamiltonians.h5"
/
get_avail_filename(
    ::Hamiltonian,
    ::Pair{Val{:source}, Val{:pred}},
    ::Pair{Val{:spin}, Val{SPIN}},
    ::Type{<:DeepHOperator}
) where SPIN = "hamiltonians_pred.h5"

get_avail_filename(
    ::Overlap, 
    ::Pair{Val{:source}, Val{:ref}},
    ::Type{<:DeepHOperator}
) = "overlaps.h5"
