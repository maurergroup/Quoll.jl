using SafeTestsets

@safetestset "common" begin

    @safetestset "atomtools.jl" begin
        include("common/atomtools.jl")
    end

    @safetestset "mpitools.jl" begin
        include("common/mpitools.jl")
    end

end

@safetestset "basis.jl" begin
    include("basis.jl")
end

@safetestset "shconversion.jl" begin
    include("shconversion.jl")
end

@safetestset "sparsity.jl" begin
    include("sparsity.jl")
end

@safetestset "spin.jl" begin
    include("spin.jl")
end

@safetestset "operators" begin

    @safetestset "canonical" begin
    
        @safetestset "bsparse.jl" begin
            include("operators/canonical/bsparse.jl")
        end

    end

    @safetestset "fhiaims" begin
    
        @safetestset "abstract.jl" begin
            include("operators/fhiaims/abstract.jl")
        end

        @safetestset "fhiaims_csc.jl" begin
            include("operators/fhiaims/fhiaims_csc.jl")
        end

    end

end

@safetestset "Parser.jl" begin
    include("Parser.jl")
end