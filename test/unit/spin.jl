using Quoll
using Dictionaries
using AtomsBase

using Test

struct DummyFormat <: Quoll.AbstractOperator end

# TODO: at the time of writing this these methods are not even in spin.jl
# so I should consider moving them out where they belong

@testset "SpinMetadata" begin
    ref = SpinMetadata(1//2)
    @test SpinMetadata(0.5) == ref

    ref = SpinMetadata(-1//2)
    @test SpinMetadata(-0.5) == ref

    @test_throws ArgumentError SpinMetadata(1)
end

@testset "SpinSetMetadata" begin
    basis = dictionary([
        ChemicalSpecies(:C) => [10, 20, 30, 40],
        ChemicalSpecies(:Si) => [10, 20, 30, 40, 50, 60],
    ])

    @testset "Spin none" begin
        spinset = Quoll.SpinSetMetadata(Hamiltonian(:foo, :none), Val(:none), basis, DummyFormat)
        @test spinset === nothing

        spinset = Quoll.SpinSetMetadata(Overlap(:foo), basis, DummyFormat)
        @test spinset === nothing
    end

    @testset "Spin up" begin
        spinset = Quoll.SpinSetMetadata(Hamiltonian(:foo, :up), Val(:up), basis, DummyFormat)
        @test spinset.soc == false
        @test spinset.spins == dictionary([
            ChemicalSpecies(:C) => fill(SpinMetadata(1//2), 4),
            ChemicalSpecies(:Si) => fill(SpinMetadata(1//2), 6),
        ])
    end
    
    @testset "Spin down" begin
        spinset = Quoll.SpinSetMetadata(Hamiltonian(:foo, :down), Val(:down), basis, DummyFormat)
        @test spinset.soc == false
        @test spinset.spins == dictionary([
            ChemicalSpecies(:C) => fill(SpinMetadata(-1//2), 4),
            ChemicalSpecies(:Si) => fill(SpinMetadata(-1//2), 6),
        ])
    end

end