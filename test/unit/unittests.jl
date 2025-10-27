using SafeTestsets

@safetestset "common" begin

    @safetestset "atomtools.jl" begin
        include("common/atomtools.jl")
    end

    @safetestset "mpitools.jl" begin
        include("common/mpitools.jl")
    end

end

@safetestset "operators" begin

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