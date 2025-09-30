module Smearing
export SmearingFunction, FermiDirac, Gaussian, implemented_smearing

abstract type SmearingFunction end

struct FermiDirac <: SmearingFunction end

struct Gaussian <: SmearingFunction end

function implemented_smearing end

implemented_smearing(::Val{:fermidirac}) = FermiDirac
implemented_smearing(::Val{:gaussian}) = Gaussian

end