using Quoll
using Configurations

using Test
using Main.TestUtils

@testset "find_leafdirs" begin

    @testset "Nominal case" begin
        tempdirs = [
            joinpath("dir1", "dir11"),
            joinpath("dir1", "dir12"),
            joinpath("dir2", "dir21"),
        ]

        root, paths = create_temptree(tempdirs)
        # root
        # ├── dir1
        # │   ├── dir11
        # │   └── dir12
        # └── dir2
        #     └── dir21

        leafdirs = Quoll.Parser.find_leafdirs(root)
        @test Set(leafdirs) == Set(paths)
    end

    @testset "Only root" begin
        tempdirs = []

        root, _ = create_temptree(tempdirs)
        # root
        # ├──

        leafdirs = Quoll.Parser.find_leafdirs(root)
        @test leafdirs == [root]

    end

end

@testset "normalize_comparison" begin

    s = "FERMI_DIRAC"
    @test Quoll.Parser.normalize_comparison(s) == "fermidirac"

    replace_pairs = tuple()
    @test Quoll.Parser.normalize_comparison(s; replace_pairs) == "fermi_dirac"

end