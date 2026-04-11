using SafeTestsets

@safetestset "tools" begin

    @safetestset "atomtools.jl" begin
        include("tools/atomtools.jl")
    end

    @safetestset "mpitools.jl" begin
        include("tools/mpitools.jl")
    end

    @safetestset "fastlookup.jl" begin
        include("tools/fastlookup.jl")
    end

    @safetestset "symmetry.jl" begin
        include("tools/symmetry.jl")
    end

end

@safetestset "core" begin

    @safetestset "components" begin

        @safetestset "basis.jl" begin
            include("core/components/basis.jl")
        end

        @safetestset "kpoints.jl" begin
            include("core/components/kpoints.jl")
        end

        @safetestset "operatorkind.jl" begin
            include("core/components/operatorkind.jl")
        end

        @safetestset "shconvention.jl" begin
            include("core/components/shconvention.jl")
        end

        @safetestset "source.jl" begin
            include("core/components/source.jl")
        end
        
        @safetestset "sparsity.jl" begin
            include("core/components/sparsity.jl")
        end

        @safetestset "spin.jl" begin
            include("core/components/spin.jl")
        end

    end
    
end

@safetestset "operators" begin
    
    @safetestset "canonical.jl" begin
        include("operators/canonical.jl")
    end
    
    @safetestset "canonical.jl" begin
        include("operators/fhiaims.jl")
    end

    @safetestset "deeph.jl" begin
        include("operators/deeph.jl")
    end

end

@safetestset "Parser.jl" begin
    include("Parser.jl")
end