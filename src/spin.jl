using Base

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

struct SpinSetMetadata
    spins::Dictionary{ChemicalSpecies, Vector{SpinMetadata}}
    soc::Bool
end
