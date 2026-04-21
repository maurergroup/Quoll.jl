using Quoll
using Documenter

DocMeta.setdocmeta!(Quoll, :DocTestSetup, :(using Quoll); recursive = true)

makedocs(;
    modules = [Quoll],
    authors = "Valdas Vitartas <valdas.vitartas@warwick.ac.uk> and contributors",
    sitename = "Quoll.jl",
    format = Documenter.HTML(;
        canonical = "https://maurergroup.github.io/Quoll.jl",
        edit_link = "main",
        assets = String[],
    ),
    pages = [
        "Home" => "index.md",
        "Getting started" => [
            "Installation" => "getting_started/installation.md",
            "What the package does" => "getting_started/what_it_does.md",
            "First example" => "getting_started/example.md",
            "Defining new methods" => "getting_started/defining_new_methods.md",
        ],
        "Tutorials" => [
            "MACE-H / DeepH-E3 training data" => "tutorials/mace_h_deeph_e3.md",
        ],
        "Input file" => "input_file.md",
        "API" => [
            "Quoll" => [
                "Core" => "api/core.md",
                "Components" => "api/components.md",
                "Operators" => "api/operators.md",
                "Tools" => "api/tools.md",
            ],
            "Quoll.Projections" => "api/projections.md",
            "Quoll.Parser" => "api/parser.md",
            "Operator interface" => "api/operator_interface.md",
        ],
        "Developer docs" => "developer.md",
    ],
)

if get(ENV, "CI", nothing) == "true"
    deploydocs(; repo = "github.com/maurergroup/Quoll.jl", devbranch = "main")
end
