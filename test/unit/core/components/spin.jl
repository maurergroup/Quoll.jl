using Quoll
using Dictionaries
using AtomsBase

using Test

# struct DummyFormat{O, T, D, M} <: Quoll.AbstractOperator{O, T, D, M} end

@testset "SpinsMetadata" begin
    a = ChemicalSpecies(:C)
    b = ChemicalSpecies(:Si)
    basis = Base.ImmutableDict(
        a => [
            Quoll.BasisMetadata(a, 1, 0, 0, nothing),
            Quoll.BasisMetadata(a, 2, 1,-1, nothing),
            Quoll.BasisMetadata(a, 2, 1, 0, nothing),
            Quoll.BasisMetadata(a, 2, 1, 1, nothing),
        ],
        b => [
            Quoll.BasisMetadata(b, 1, 0, 0, nothing),
            Quoll.BasisMetadata(b, 2, 0, 0, nothing),
            Quoll.BasisMetadata(b, 3, 0, 0, nothing),
            Quoll.BasisMetadata(b, 2, 1,-1, nothing),
            Quoll.BasisMetadata(b, 2, 1, 0, nothing),
            Quoll.BasisMetadata(b, 2, 1, 1, nothing),
        ],
    )

    # @testset "Spin none" begin
    #     spins = Quoll.SpinsMetadata(Hamiltonian(source = :foo, spin = :none), Val(:none), basis, DummyFormat)
    #     @test spins === nothing

    #     spins = Quoll.SpinsMetadata(Overlap(source = :foo), basis, DummyFormat)
    #     @test spins === nothing
    # end

    # @testset "Spin up" begin
    #     spins = Quoll.SpinsMetadata(Hamiltonian(source = :foo, spin = :up), Val(:up), basis, DummyFormat)
    #     @test spins.soc == false
    #     @test spins.spins == Base.ImmutableDict(
    #         ChemicalSpecies(:C) => fill(⬆, 4),
    #         ChemicalSpecies(:Si) => fill(⬆, 6),
    #     )
    # end
    
    # @testset "Spin down" begin
    #     spins = Quoll.SpinsMetadata(Hamiltonian(source = :foo, spin = :down), Val(:down), basis, DummyFormat)
    #     @test spins.soc == false
    #     @test spins.spins == Base.ImmutableDict(
    #         ChemicalSpecies(:C) => fill(⬇, 4),
    #         ChemicalSpecies(:Si) => fill(⬇, 6),
    #     )
    # end

end