abstract type AbstractFHIaimsOperator{O, T, D, M} <: AbstractOperator{O, T, D, M} end
abstract type AbstractFHIaimsMetadata{A, S, B, P} <: AbstractOperatorMetadata{A, S, B, P} end

# Same order as wiki
# For m > 0 opposite CS phase
const FHIaimsSHConversion = SHConversion(
    [[1], [1,  2,  3], [1,  2,  3,  4,  5], [1,  2,  3,  4,  5,  6,  7], [1,  2,  3,  4,  5,  6,  7,  8,  9]],
    [[1], [1,  1, -1], [1,  1,  1, -1,  1], [1,  1,  1,  1, -1,  1, -1], [1,  1,  1,  1,  1, -1,  1, -1,  1]],
)

SHConversion(::Type{<:AbstractFHIaimsOperator}) = FHIaimsSHConversion

function load_atoms(dir::AbstractString, ::Type{<:AbstractFHIaimsOperator})
    path = joinpath(dir, "geometry.in")
    atoms = load_system(path)

    # Recentre atom positions to be consistent with relative positions
    # that are used during real-space matrix integration in FHI-aims
    all(periodicity(atoms)) && (atoms = recentre(atoms))
    
    return atoms
end

function BasisSetMetadata(dir::AbstractString, atoms::AbstractSystem, ::Type{<:AbstractFHIaimsOperator})
    p = joinpath(dir, "basis-indices.out")
    @argcheck ispath(p)

    atom2species = species(atoms, :)
    basis = Base.ImmutableDict((
        species => BasisMetadata{Base.ImmutableDict{String, String}}[]
        for species in unique(atom2species)
    )...)

    f = open(p, "r")

    # Reads both the empty lines and the header
    while isempty(strip(readline(f))) end

    species2firstatom = dictionary(
        species => findfirst(x -> x == species, atom2species)
        for species in unique(atom2species)
    )

    ib = 1
    while !(eof(f))
        _, type, iat, n, l, m = convert.(String, split(readline(f)))
        iat, n, l, m = parse.(Int, (iat, n, l, m))

        if iat == species2firstatom[atom2species[iat]]
            push!(
                basis[atom2species[iat]],
                BasisMetadata(atom2species[iat], n, l, m, Base.ImmutableDict("type" => type))
            )
        end

        ib += 1
    end
    close(f)

    return BasisSetMetadata(basis, atom2species)
end