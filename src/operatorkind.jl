abstract type AbstractOperatorKind end

struct Hamiltonian <:AbstractOperatorKind
    tag::Symbol
end

struct Overlap <:AbstractOperatorKind
    tag::Symbol
end

function get_operatorkind end

get_operatorkind(::Val{:h}) = Hamiltonian(:ref)
get_operatorkind(::Val{:href}) = Hamiltonian(:ref)
get_operatorkind(::Val{:hpred}) = Hamiltonian(:pred)

get_operatorkind(::Val{:s}) = Overlap(:ref)
get_operatorkind(::Val{:sref}) = Overlap(:ref)
get_operatorkind(::Val{:spred}) = Overlap(:pred)
