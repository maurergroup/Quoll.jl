using Quoll
using Configurations
using Test

include("../testutils.jl")

@testset "find_leafdirs" begin

    tempdirs = [
        joinpath("dir", "dir1", "dir11"),
        joinpath("dir", "dir1", "dir12"),
        joinpath("dir", "dir2", "dir21"),
    ]

    root, paths = create_tempdirs(tempdirs)
    # dir
    # ├── dir1
    # │   ├── dir11
    # │   └── dir12
    # └── dir2
    #     └── dir21

    leafdirs = Quoll.Utils.find_leafdirs(root)
    @test Set(leafdirs) == Set(paths)
end

@testset "normalize_comparison" begin

    s = "FERMI_DIRAC"
    @test Quoll.Utils.normalize_comparison(s) == "fermidirac"

    replace_pairs = tuple()
    @test Quoll.Utils.normalize_comparison(s; replace_pairs) == "fermi_dirac"

end