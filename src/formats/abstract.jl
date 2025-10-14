abstract type AbstractOperator end
abstract type AbstractOperatorMetadata end

# TODO: Even though the fields of metadata might be the same among different formats,
# the constructors will differ.
# An alternative would be to use dispatch on the format instead,
# e.g. OperatorMetadata(..., ::Type{<:AbstractOperator}).
# For now I will keep metadata as separate structs for each format,
# but if I find that fields are the same I should consider refactoring
# the code in favour of dispatch on operator format (more consistent
# with the rest and fewer unnecessary types).

function get_readformat end
function get_writeformat end

function get_avail_operatorkinds end
function get_avail_filenames end

# Find available operators in the supplied dir
function find_operatorkinds(dir::AbstractString, format::Type{<:AbstractOperator})
    found_operatorkinds = AbstractOperatorKind[]
    for operatorkind in get_avail_operatorkinds(format)
        names = get_avail_filenames(operatorkind, format)
        ps = joinpath.(dir, names)
        if any(ispath.(ps))
            push!(found_operatorkinds, operatorkind)
        end
    end
    return found_operatorkinds
end

function get_avail_filenames(operatorkind::Hamiltonian, T::Type{<:AbstractOperator})
    return get_avail_filenames(operatorkind, Val(operatorkind.source), Val(operatorkind.spin), T)
end

function get_avail_filenames(operatorkind::Overlap, T::Type{<:AbstractOperator})
    return get_avail_filenames(operatorkind, Val(operatorkind.source), T)
end

# Constructors for SpinSetMetadata from AbstractOperatorKind, BasisSetMetadata, and AbstractOperator.
# AbstractOperator is required for the SOC case. This is because for now in the case of SOC
# BasisSetMetadata.basis stores 2x the number of basis functions, and we don't know which of them
# are spin up and spin down. This type of information might depend on the format
# TODO: this could be moved to spin.jl if this file was imported before spin.jl
# (the only reason this is here is because of the dependence on AbstractOperator)

function SpinSetMetadata(kind::AbstractOperatorKind, basisset::BasisSetMetadata, format::Type{<:AbstractOperator})
    return SpinSetMetadata(kind, basisset.basis, format)
end

function SpinSetMetadata(kind::Hamiltonian, basis::Dictionary{ChemicalSpecies}, format::Type{<:AbstractOperator})
    return SpinSetMetadata(kind, Val(kind.spin), basis, format)
end

function SpinSetMetadata(::Hamiltonian, ::Val{:up}, basis::Dictionary{ChemicalSpecies}, ::Type{<:AbstractOperator})
    return SpinSetMetadata(
        Dictionary(
            keys(basis),
            [fill(SpinMetadata(1//2), length(atom_basis)) for atom_basis in basis]
        ),
        false,
    )
end

function SpinSetMetadata(::Hamiltonian, ::Val{:down}, basis::Dictionary{ChemicalSpecies}, ::Type{<:AbstractOperator})
    return SpinSetMetadata(
        Dictionary(
            keys(basis),
            [fill(SpinMetadata(-1//2), length(atom_basis)) for atom_basis in basis]
        ),
        false,
    )
end

SpinSetMetadata(::Hamiltonian, ::Val{:none}, ::Dictionary{ChemicalSpecies}, ::Type{<:AbstractOperator}) = nothing

SpinSetMetadata(::Overlap, ::Dictionary{ChemicalSpecies}, ::Type{<:AbstractOperator}) = nothing