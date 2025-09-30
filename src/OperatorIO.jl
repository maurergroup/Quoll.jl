module OperatorIO
export implemented_read_format, implemented_write_format
export AtomBasisOperator
export FHIaimsOperator, DeepHOperator

abstract type AtomBasisOperator end

abstract type AbstractQuollOperator <: AtomBasisOperator end
abstract type AbstractBSparseOperator <: AbstractQuollOperator end
struct RealBSparseOperator <: AbstractBSparseOperator end

# This will not be a singleton, which means if we want to dispatch
# on the format without passing values into a method we will need to dispatch
# on the type instead of an instance (::Type{FHIaimsOperator}).
# Would that be inconsistent if I use singleton instances for e.g. FC99V?
struct FHIaimsOperator <: AtomBasisOperator end

struct DeepHOperator <: AtomBasisOperator end

function implemented_read_format end
function implemented_write_format end

# Implemented read formats
implemented_read_format(::Val{:fhiaims}) = FHIaimsOperator

# Implemented write formats
implemented_write_format(::Val{:deeph}) = DeepHOperator

end