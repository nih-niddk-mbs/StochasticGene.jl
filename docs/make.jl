push!(LOAD_PATH,"../src/")
using StochasticGene
using Documenter
makedocs(
         sitename = "StochasticGene.jl",
         modules  = [StochasticGene],
         pages=[
                "Home" => "index.md"
               ])
deploydocs(;
    repo="github.com/nih-niddk-mbs/StochasticGene.jl.git",
)
