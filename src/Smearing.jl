module Smearing
export SmearingFunction, FermiDirac, Gaussian, SMEAR_REGISTRY

abstract type SmearingFunction end

struct FermiDirac <: SmearingFunction end

struct Gaussian <: SmearingFunction end

const SMEAR_REGISTRY = (FermiDirac, Gaussian)

end