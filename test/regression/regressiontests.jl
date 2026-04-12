using SafeTestsets

@safetestset "Silicon carbide FHIaimsCSC -> CanonicalBlock -> DeepHBlock" begin
    include("fhi_canon_deeph_siliconcarbide.jl")
end

# @safetestset "Gold FHIaimsCSC -> CanonicalBlock -Laikov-> CanonicalBlock -> DeepHBlock" begin
#     include("fhi_canon_laikov_deeph_gold.jl")
# end

# @safetestset "N-doped graphene DeepHBlock -> CanonicalBlock -> CanonicalDenseRecip" begin
#     include("deeph_canon_canonrecip_dopedgraphene.jl")
# end
