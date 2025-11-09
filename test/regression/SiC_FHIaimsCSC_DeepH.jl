using Test
using Main.TestUtils

np = 2
project = joinpath(@__DIR__, "../../app")
app = joinpath(@__DIR__, "../../app/quoll.jl")

inputfile = """
[input]
format = "FHI-aims"
directory = "SiC_FHIaims/SiC_FHIaims"
operators = ["H", "S"]

[output]
format = "DeepH"
directory = "SiC_DeepH_converted"
"""

SiC_examples = joinpath(@__DIR__, "../../examples/SiC")
tarballs = [
    joinpath(SiC_examples, "SiC_FHIaims.tar.gz")
    joinpath(SiC_examples, "SiC_DeepH.tar.gz")
]

setupteardown_tmp(tarballs = tarballs) do
    inputfile_path = "inputfile.dat"
    write(inputfile_path, inputfile)
    run(mpiexec_quollapp(app, project, inputfile_path; np = np))
    @test true
end