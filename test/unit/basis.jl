using Quoll
using Dictionaries
using AtomsBase

using Test

const basisset = BasisSetMetadata(
    dictionary([
        ChemicalSpecies(:H1) => [
            Quoll.BasisMetadata(a, 1, 0, 0, nothing)
            Quoll.BasisMetadata(a, 2, 0, 0, nothing)
            Quoll.BasisMetadata(a, 2, 1,-1, nothing)
            Quoll.BasisMetadata(a, 2, 1, 0, nothing)
            Quoll.BasisMetadata(a, 2, 1, 1, nothing)
        ],
        ChemicalSpecies(:H2) => [
            Quoll.BasisMetadata(b, 2, 0, 0, nothing)
            Quoll.BasisMetadata(b, 2, 1,-1, nothing)
            Quoll.BasisMetadata(b, 2, 1, 0, nothing)
            Quoll.BasisMetadata(b, 2, 1, 1, nothing)
        ]
    ]),
    [ChemicalSpecies(:H1), ChemicalSpecies(:H2), ChemicalSpecies(:H2), ChemicalSpecies(:H1)]
)

@testset "get_atom2nbasis" begin
    ref_atom2nbasis = [5, 4, 4, 5]
    @test ref_atom2nbasis == get_atom2nbasis(basisset)
end

@testset "get_atom2basis" begin
    ref_atom2basis = [1:5, 6:9, 10:14, 15:18]
    @test ref_atom2basis == get_atom2basis(basisset)
end

@testset "get_basis2atom" begin
    ref_basis2atom = [1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4]
    @test ref_basis2atom == get_basis2atom(basisset)
end

@testset "get_species2nbasis" begin
    a = ChemicalSpecies(:H1)
    b = ChemicalSpecies(:H2)

    ref_species2nbasis = dictionary([a => 5, b => 4])
    species2nbasis = Quoll.get_species2nbasis(basisset)

    @test keys(species2nbasis) == keys(ref_species2nbasis)
    for z in keys(ref_species2nbasis)
        @test species2nbasis[z] == ref_species2nbasis[z]
    end
end

@testset "basis_atom" begin
    a = ChemicalSpecies(:H1)
    b = ChemicalSpecies(:H2)
    basisset = basisset

    basis_atom3 = [
        Quoll.BasisMetadata(b, 2, 0, 0, nothing)
        Quoll.BasisMetadata(b, 2, 1,-1, nothing)
        Quoll.BasisMetadata(b, 2, 1, 0, nothing)
        Quoll.BasisMetadata(b, 2, 1, 1, nothing)
    ]
    @test Quoll.basis_atom(basisset, 3) == basis_atom3

    basis_atom4 = [
        Quoll.BasisMetadata(a, 1, 0, 0, nothing)
        Quoll.BasisMetadata(a, 2, 0, 0, nothing)
        Quoll.BasisMetadata(a, 2, 1,-1, nothing)
        Quoll.BasisMetadata(a, 2, 1, 0, nothing)
        Quoll.BasisMetadata(a, 2, 1, 1, nothing)
    ]
    @test Quoll.basis_atom(basisset, 4) == basis_atom4
end

@testset "basis_species" begin
    a = ChemicalSpecies(:H1)
    b = ChemicalSpecies(:H2)
    basisset = basisset

    basis_species_b = [
        Quoll.BasisMetadata(b, 2, 0, 0, nothing)
        Quoll.BasisMetadata(b, 2, 1,-1, nothing)
        Quoll.BasisMetadata(b, 2, 1, 0, nothing)
        Quoll.BasisMetadata(b, 2, 1, 1, nothing)
    ]
    @test Quoll.basis_species(basisset, b) == basis_species_b

    basis_species_a = [
        Quoll.BasisMetadata(a, 1, 0, 0, nothing)
        Quoll.BasisMetadata(a, 2, 0, 0, nothing)
        Quoll.BasisMetadata(a, 2, 1,-1, nothing)
        Quoll.BasisMetadata(a, 2, 1, 0, nothing)
        Quoll.BasisMetadata(a, 2, 1, 1, nothing)
    ]
    @test Quoll.basis_species(basisset, a) == basis_species_a
end

@testset "get_offsets" begin
    ref_offsets = [0, 5, 9, 14]
    @test ref_offsets == Quoll.get_offsets(basisset)
end