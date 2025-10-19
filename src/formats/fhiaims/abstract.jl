abstract type AbstractFHIaimsOperator <: AbstractOperator end
abstract type AbstractFHIaimsMetadata <: AbstractOperatorMetadata end

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

    n_atoms = length(atoms)
    n_basis_atom = zeros(Int, n_atoms)
    atom2species = species(atoms, :)
    basis2atom = Int[]

    basis = dictionary(
        species => BasisMetadata{Dict{String, String}}[]
        for species in unique(atom2species)
    )

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
                BasisMetadata(atom2species[iat], n, l, m, Dict{String, String}("type" => type))
            )
        end

        push!(basis2atom, iat)
        n_basis_atom[iat] += 1
        ib += 1
    end
    close(f)

    atom2basis_start = cumsum(n_basis_atom) .- n_basis_atom .+ 1
    atom2basis = [
        atom2basis_start[iat]:atom2basis_start[iat] + n_basis_atom[iat] - 1
        for iat in 1:n_atoms
    ]

    return BasisSetMetadata(basis, n_basis_atom, atom2species, basis2atom, atom2basis)
end