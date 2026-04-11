using Quoll

using Test

struct DummySource <: Quoll.AbstractSource
    extra_param::Symbol
end

@testset "namedtuple" begin
    source = DummySource(:a)
    namedtuple_ref = (; extra_param = :a)
    @test Quoll.namedtuple(source) == namedtuple_ref
end