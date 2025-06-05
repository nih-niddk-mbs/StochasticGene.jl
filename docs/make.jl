using Documenter
using StochasticGene
using Pkg

# Read version from main Project.toml for synchronization
project_path = abspath(joinpath(@__DIR__, "..", "Project.toml"))
pkg_version = Pkg.TOML.parsefile(project_path)["version"]

makedocs(
    sitename = "StochasticGene.jl",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://nih-niddk-mbs.github.io/StochasticGene.jl/stable/",
        assets = ["assets/custom.css"],
        analytics = "UA-XXXXXXXXX-X",
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
    version = pkg_version,  # Always synchronized with package
    doctest = false,        # Disable running doctests
    linkcheck = false,      # Disable link checking
    checkdocs = :none,      # Disable documentation checking
    warnonly = true         # Only show warnings, don't fail on warnings
)

deploydocs(
    repo = "github.com/nih-niddk-mbs/StochasticGene.jl.git",
    devbranch = "main",
    target = "build",
    branch = "gh-pages",
    push_preview = true,
    forcepush = true,
    versions = ["stable" => "v^", "v#.#", "dev" => "dev"],
    make = nothing,
)
