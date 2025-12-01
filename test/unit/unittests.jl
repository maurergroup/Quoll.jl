using SafeTestsets

@safetestset "methods" begin

    @safetestset "atomtools.jl" begin
        include("methods/atomtools.jl")
    end

    @safetestset "mpitools.jl" begin
        include("methods/mpitools.jl")
    end

    @safetestset "fastlookup.jl" begin
        include("methods/fastlookup.jl")
    end

    @safetestset "symmetry.jl" begin
        include("methods/symmetry.jl")
    end

end

@safetestset "basis.jl" begin
    include("basis.jl")
end

@safetestset "kpoints.jl" begin
    include("kpoints.jl")
end

@safetestset "operatorkind.jl" begin
    include("operatorkind.jl")
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