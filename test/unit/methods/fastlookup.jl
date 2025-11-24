using Quoll
using AtomsBase

using Test

const a = ChemicalSpecies(:H1)
const b = ChemicalSpecies(:H2)

const atom2species = begin
    [a, b, b]
end

@testset "SpeciesDict{Vector{Int}}" begin
    d = Dict(
        a => [10, 20],
        b => [30, 40]
    )

    ref_atomarray = Vector{Vector{Int}}(undef, 3)
    ref_atomarray[1] = [10, 20]
    ref_atomarray[2] = [30, 40]
    ref_atomarray[3] = [30, 40]

    atomarray = Quoll.convert_to_atomarray(d, atom2species)
    @test atomarray == ref_atomarray
end

@testset "SpeciesPairDict{Vector{Int}}" begin
    d = Dict(
        (a, a) => [10, 20],
        (a, b) => [30, 40],
        (b, b) => [50, 60],
    )

    ref_atomarray = Matrix{Vector{Int}}(undef, 3, 3)
    ref_atomarray[1, 1] = Int[10, 20]
    ref_atomarray[1, 2] = Int[30, 40]
    ref_atomarray[1, 3] = Int[30, 40]
    ref_atomarray[2, 1] = Int[]
    ref_atomarray[2, 2] = Int[50, 60]
    ref_atomarray[2, 3] = Int[50, 60]
    ref_atomarray[3, 1] = Int[]
    ref_atomarray[3, 2] = Int[50, 60]
    ref_atomarray[3, 3] = Int[50, 60]

    atomarray = Quoll.convert_to_atomarray(d, atom2species)
    @test atomarray == ref_atomarray
end

@testset "AtomPairAnyDict{ThreeDimSliceView{Float64}}" begin
    sliceview = view(fill(1.0, 2, 2, 4), :, :, 1:2)
    sentinel = view(zeros(Float64, 1, 1, 1), :, :, 0:-1)

    d = Dict(
        (1, 1) => sliceview,
        (1, 2) => sliceview,
        (3, 3) => sliceview
    )

    ref_atomarray = Matrix{typeof(sliceview)}(undef, 3, 3)
    ref_atomarray[1, 1] = sliceview
    ref_atomarray[1, 2] = sliceview
    ref_atomarray[1, 3] = sentinel
    ref_atomarray[2, 1] = sentinel
    ref_atomarray[2, 2] = sentinel
    ref_atomarray[2, 3] = sentinel
    ref_atomarray[3, 1] = sentinel
    ref_atomarray[3, 2] = sentinel
    ref_atomarray[3, 3] = sliceview

    atomarray = Quoll.convert_to_atomarray(d, 3)
    @test atomarray == ref_atomarray
end