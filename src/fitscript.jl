# @everywhere include("StochasticGene.jl")
#
# using Dates

include("/home/carsonc/StochasticGene/src/fit_rna_T120.jl")

gene = ARGS[1]

fit_rna(gene)
