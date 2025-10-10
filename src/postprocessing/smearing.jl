abstract type SmearingFunction end

struct FermiDirac <: SmearingFunction end
struct Gaussian <: SmearingFunction end

function get_smearing end

get_smearing(::Val{:fermidirac}) = FermiDirac
get_smearing(::Val{:gaussian}) = Gaussian
