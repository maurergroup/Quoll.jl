using SafeTestsets

@safetestset "Silicon carbide FHIaimsCSC -> CanonicalBlock -> DeepHBlock" begin
    include("fhi_canon_deeph_siliconcarbide.jl")
end

@safetestset "Gold FHIaimsCSC -> CanonicalBlock -Laikov-> CanonicalBlock -> DeepHBlock" begin
    include("fhi_canon_laikov_deeph_gold.jl")
end

@safetestset "Water FHIaimsCSC -> CanonicalBlock -Laikov-> CanonicalBlock -> DeepHBlock" begin
    include("fhi_canon_laikov_deeph_water.jl")
end

@safetestset "Bismuth telluride DeepHBlockSpin -> CanonicalBlockSpin -fourier,inv_fourier-> CanonicalBlockSpin -> DeepHBlockSpin" begin
    include("deeph_canon_recip_bismuthtelluride.jl")
end

@safetestset "Graphene DeepHBlock,FHIaimsCSC -> CanonicalBlock -zero_out-> CanonicalBlock" begin
    include("deeph_canon_zeroout_graphene.jl")
end
