abstract type AbstractOperatorKind end

struct Hamiltonian <:AbstractOperatorKind
    tag::Symbol
end

struct Overlap <:AbstractOperatorKind
    tag::Symbol
end

operator_kind(::Val{:h}) = Hamiltonian(:ref)
operator_kind(::Val{:href}) = Hamiltonian(:ref)
operator_kind(::Val{:hpred}) = Hamiltonian(:pred)

operator_kind(::Val{:s}) = Overlap(:ref)
operator_kind(::Val{:sref}) = Overlap(:ref)
operator_kind(::Val{:spred}) = Overlap(:pred)
