using Quoll
using Documenter

DocMeta.setdocmeta!(Quoll, :DocTestSetup, :(using Quoll); recursive=true)

makedocs(;
    modules=[Quoll],
    authors="Valdas Vitartas <valdas.vitartas@warwick.ac.uk> and contributors",
    sitename="Quoll.jl",
    format=Documenter.HTML(;
        canonical="https://maurergroup.github.io/Quoll.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/maurergroup/Quoll.jl",
    devbranch="main",
)
