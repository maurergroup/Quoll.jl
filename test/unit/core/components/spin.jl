using Quoll
using Dictionaries
using AtomsBase

using Test

struct DummySource <: Quoll.AbstractSource end

const basisset = let a = ChemicalSpecies(:H1), b = ChemicalSpecies(:H2)
    Quoll.BasisSetMetadata(
        dictionary([
            a => [
                Quoll.BasisMetadata(a, 1, 0, 0, nothing)
                Quoll.BasisMetadata(a, 2, 0, 0, nothing)
                Quoll.BasisMetadata(a, 2, 1, -1, nothing)
                Quoll.BasisMetadata(a, 2, 1, 0, nothing)
                Quoll.BasisMetadata(a, 2, 1, 1, nothing)
            ],
            b => [
                Quoll.BasisMetadata(b, 2, 0, 0, nothing)
                Quoll.BasisMetadata(b, 2, 1, -1, nothing)
                Quoll.BasisMetadata(b, 2, 1, 0, nothing)
                Quoll.BasisMetadata(b, 2, 1, 1, nothing)
            ],
        ]),
        [a, b, b, a],
    )
end

const basisset_soc = let a = ChemicalSpecies(:H1), b = ChemicalSpecies(:H2)
    Quoll.BasisSetMetadata(
        dictionary([
            a => [
                Quoll.BasisMetadata(a, 1, 0, 0, nothing)
                Quoll.BasisMetadata(a, 2, 0, 0, nothing)
                Quoll.BasisMetadata(a, 2, 1, -1, nothing)
                Quoll.BasisMetadata(a, 2, 1, 0, nothing)
                Quoll.BasisMetadata(a, 2, 1, 1, nothing)
                Quoll.BasisMetadata(a, 1, 0, 0, nothing)
                Quoll.BasisMetadata(a, 2, 0, 0, nothing)
                Quoll.BasisMetadata(a, 2, 1, -1, nothing)
                Quoll.BasisMetadata(a, 2, 1, 0, nothing)
                Quoll.BasisMetadata(a, 2, 1, 1, nothing)
            ],
            b => [
                Quoll.BasisMetadata(b, 2, 0, 0, nothing)
                Quoll.BasisMetadata(b, 2, 1, -1, nothing)
                Quoll.BasisMetadata(b, 2, 1, 0, nothing)
                Quoll.BasisMetadata(b, 2, 1, 1, nothing)
                Quoll.BasisMetadata(b, 2, 0, 0, nothing)
                Quoll.BasisMetadata(b, 2, 1, -1, nothing)
                Quoll.BasisMetadata(b, 2, 1, 0, nothing)
                Quoll.BasisMetadata(b, 2, 1, 1, nothing)
            ],
        ]),
        [a, b, b, a],
    )
end

@testset "SpinsMetadata" begin
    a = ChemicalSpecies(:H1)
    b = ChemicalSpecies(:H2)

    @testset "Spin up" begin
        spins = Quoll.SpinsMetadata(
            DummySource(),
            Hamiltonian(source = :foo, spin = :up),
            basisset,
        )
        @test spins.soc == false
        @test spins.spins == Base.ImmutableDict(
            a => fill(Quoll.⬆, 5),
            b => fill(Quoll.⬆, 4),
        )
    end

    @testset "Spin down" begin
        spins = Quoll.SpinsMetadata(
            DummySource(),
            Hamiltonian(source = :foo, spin = :down),
            basisset,
        )
        @test spins.soc == false
        @test spins.spins == Base.ImmutableDict(
            a => fill(Quoll.⬇, 5),
            b => fill(Quoll.⬇, 4),
        )
    end

    @testset "SOC" begin
        spins = Quoll.SpinsMetadata(
            DummySource(),
            Hamiltonian(source = :foo, spin = :soc),
            basisset_soc,
        )
        @test spins.soc == true
        @test spins.spins == Base.ImmutableDict(
            a => vcat(fill(Quoll.⬆, 5), fill(Quoll.⬇, 5)),
            b => vcat(fill(Quoll.⬆, 4), fill(Quoll.⬇, 4)),
        )
    end
end

@testset "reduce_spins" begin
    a = ChemicalSpecies(:H1)
    b = ChemicalSpecies(:H2)
    subbasis = [
        Quoll.BasisMetadata(a, 2, 1, -1, nothing),
        Quoll.BasisMetadata(a, 2, 1, 0, nothing),
        Quoll.BasisMetadata(a, 2, 1, 1, nothing),
        Quoll.BasisMetadata(b, 2, 0, 0, nothing)
    ]

    spins = Quoll.SpinsMetadata(
        DummySource(),
        Quoll.Hamiltonian(source = :ref, spin = :up),
        basisset,
    )
    spins_reduced = Quoll.reduce_spins(spins, basisset, subbasis)
    spins_reduced_ref = Quoll.SpinsMetadata(
        Base.ImmutableDict(
            a => fill(Quoll.⬆, 3),
            b => fill(Quoll.⬆, 1),
        ),
        false,
    )
    @test spins_reduced.soc == spins_reduced_ref.soc
    @test spins_reduced.spins == spins_reduced_ref.spins

end