using SafeTestsets

# TODO: Could separate into two different tests when RealBSparse
# contains write methods
@safetestset "SiC FHIaimsCSC -> RealBSparse -> DeepH" begin
    include("SiC_FHIaimsCSC_DeepH.jl")
end
