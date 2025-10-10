struct DeepHOperator <: AbstractOperator end

get_writeformat(::Val{:deeph}) = DeepHOperator