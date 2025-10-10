abstract type SmearingFunction end

struct FermiDirac <: SmearingFunction end
struct Gaussian <: SmearingFunction end

function smearing end

smearing(::Val{:fermidirac}) = FermiDirac
smearing(::Val{:gaussian}) = Gaussian
