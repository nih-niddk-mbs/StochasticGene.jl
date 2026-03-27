@everywhere using StochasticGene
@time fit(; key="3311", datapath="data/5Prime_gene_enhancer/including_background/short", maxtime=30.0, root=".", resultfolder="coupled", infolder="coupled", basekeys="3311")
