using Documenter
using StochasticGene

makedocs(
    sitename = "StochasticGene.jl",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", "false") == "true",
        assets = ["assets/custom.css"],
    ),
    modules = [StochasticGene],
    pages = [
        "Home" => "index.md",
        "Installation" => "installation.md",
        "Getting Started" => "getting_started.md",
        "API" => [
            "API Reference" => "api/index.md"
        ],
        "Contributing" => "contributing.md",
    ],
    doctest = false,
    linkcheck = false,
    checkdocs = :none,
)

deploydocs(
    repo = "github.com/nih-niddk-mbs/StochasticGene.jl.git",
    devbranch = "main",
    push_preview = true,
    forcepush = true,
)
