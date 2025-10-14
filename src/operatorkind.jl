abstract type AbstractOperatorKind end

struct Hamiltonian <:AbstractOperatorKind
    source::Symbol
    spin::Symbol
end

Hamiltonian(source::Symbol) = Hamiltonian(source, :none)

struct Overlap <:AbstractOperatorKind
    source::Symbol
end

function get_operatorkinds end

get_operatorkinds(::Val{:h}) = [
    Hamiltonian(source, spin)
    for source in (:ref, :pred) for spin in (:none, :soc, :up, :down)
]
get_operatorkinds(::Val{:href}) = [
    Hamiltonian(source, spin)
    for source in (:ref,) for spin in (:none, :soc, :up, :down)
]
get_operatorkinds(::Val{:hpred}) = [
    Hamiltonian(source, spin)
    for source in (:pred,) for spin in (:none, :soc, :up, :down)
]
get_operatorkinds(::Val{:hsoc}) = [
    Hamiltonian(source, spin)
    for source in (:ref, :pred) for spin in (:soc,)
]
get_operatorkinds(::Val{:hpol}) = [
    Hamiltonian(source, spin)
    for source in (:ref, :pred) for spin in (:up, :down)
]

get_operatorkinds(::Val{:s}) = [
    Overlap(source)
    for source in (:ref, :pred)
]
get_operatorkinds(::Val{:sref}) = [
    Overlap(source)
    for source in (:ref,)
]
get_operatorkinds(::Val{:spred}) = [
    Overlap(source)
    for source in (:pred,)
]
