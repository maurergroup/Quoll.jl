using AutoHashEquals
using Base

@auto_hash_equals struct BasisMetadata{E}
    z::ChemicalSpecies
    n::Int
    l::Int
    m::Int
    extras::E
end
function Base.show(io::IO, basisf::BasisMetadata)
    print(io, "$(basisf.z)$(basisf.n)$(basisf.l)$(basisf.m)")
end

struct BasisSetMetadata{E}
    basis::Dictionary{ChemicalSpecies, Vector{BasisMetadata{E}}}
    n_basis_atom::Vector{Int}
    atom2species::Vector{ChemicalSpecies}
    basis2atom::Vector{Int}
    atom2basis::Vector{UnitRange{Int}}
    # TODO: spherical harmonics convention? Can be useful to have generally.
    # However, I'm not sure if I could use it to perform conversion between
    # different conventions in an efficient way
end

struct SpinMetadata
    spin::Rational{Int}

    function SpinMetadata(spin)
        @argcheck spin ∈ (-1//2, 1//2)
        return new(spin)
    end
end

function Base.show(io::IO, spin::SpinMetadata)
    if spin.spin == -1//2
        print(io, "↓")
    elseif spin.spin == 1//2
        print(io, "↑")
    end
end