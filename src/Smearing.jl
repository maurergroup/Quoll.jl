module Smearing
export SmearingFunction, FermiDirac, Gaussian, smearing

abstract type SmearingFunction end

struct FermiDirac <: SmearingFunction end

struct Gaussian <: SmearingFunction end

function smearing end

# TODO: move to parser? Not sure...
smearing(::Val{:fermidirac}) = FermiDirac
smearing(::Val{:gaussian}) = Gaussian

end