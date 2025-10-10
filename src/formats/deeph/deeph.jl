struct DeepHOperator <: AbstractOperator end

write_format(::Val{:deeph}) = DeepHOperator