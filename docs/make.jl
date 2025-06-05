using Documenter
using StochasticGene

# Generate documentation
makedocs(
    sitename = "StochasticGene.jl",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://nih-niddk-mbs.github.io/StochasticGene.jl/stable/",
        assets = String[],
        analytics = "UA-XXXXXXXXX-X",
    ),
    modules = [StochasticGene],
    pages = [
        "Home" => "index.md",
        "Installation" => "installation.md",
        "Getting Started" => "getting_started.md",
        "API" => [
            "API Reference" => "api/index.md",
        ],
        "Examples" => "examples.md",
        "Contributing" => "contributing.md",
    ],
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/nih-niddk-mbs/StochasticGene.jl.git",
    devbranch = "main",
    push_preview = true,
    versions = ["stable" => "v^", "v#.#", "dev" => "dev"],
) 