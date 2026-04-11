using Quoll
using Dictionaries

using Test

@testset "OperatorKind" begin
    h1 = Hamiltonian(source = :ref, spin = :none)
    @test h1 isa Quoll.OperatorKind{:Hamiltonian}
    @test h1.source == :ref
    @test h1.spin == :none

    h2 = Hamiltonian(Pair(:source, :ref), Pair(:spin, :none))
    @test h2.source == :ref
    @test h2.spin == :none
end

@testset "allows_symmetry" begin
    h1 = Hamiltonian(source = :ref, spin = :none)
    @test Quoll.allows_symmetry(h1) == true

    h2 = Hamiltonian(source = :ref, spin = :soc)
    @test Quoll.allows_symmetry(h2) == false

    s = Overlap(source = :ref)
    operatorkinds = [h2, s]
    @test Quoll.allows_symmetry(operatorkinds) == false
end

@testset "get_tags" begin
    h1 = Hamiltonian(source = :ref, spin = :none)
    @test Quoll.get_tags(h1, excluded_keys = nothing) == UnorderedDictionary([:source, :spin], [:ref, :none])
    @test Quoll.get_tags(h1, excluded_keys = :spin) == UnorderedDictionary([:source], [:ref])
end

@testset "get_operator_groups" begin
    h1 = Hamiltonian(source = :ref, spin = :none)
    h2 = Hamiltonian(source = :pred, spin = :none)
    h3 = Hamiltonian(source = :ref, spin = :up)
    h4 = Hamiltonian(source = :ref, spin = :down)
    s1 = Overlap(source = :ref, spin = :none)
    s2 = Overlap(source = :pred, spin = :none)
    ops = [h1, h2, h3, h4, s1, s2]

    @test Quoll.get_operator_groups(ops; excluded_keys = [:source]) == [[h1, h2, s1, s2], [h3], [h4]]
    @test Quoll.get_operator_groups(ops; op_filter = op -> op.spin == :none) == [[h1, s1], [h2, s2]]
end

@testset "sift_tags" begin
    tags = (:A, :B, :C)
    @test Quoll.sift_tags(tags) == (:A, :B, :C)

    @test Quoll.sift_tags(tags; requested_tags = (:B,)) == (:B,)
    @test_throws ErrorException Quoll.sift_tags(tags, requested_tags = (:D,))

    @test Quoll.sift_tags(tags; removed_tags = (:B,)) == (:A, :C)
    @test Quoll.sift_tags(tags; removed_tags = (:D,)) == (:A, :B, :C)
end