# Script executed by Biowulf to run code

@everywhere include("/home/carsonc/StochasticGene/src/StochasticGene.jl")

include("/home/carsonc/StochasticGene/runfiles/scriptfunctions.jl")

@time fit_rna(parse(Int,ARGS[1]),ARGS[2],ARGS[3],parse(Int,ARGS[4]),parse(Float64,ARGS[5]),ARGS[6],ARGS[7],ARGS[8],ARGS[9],ARGS[10],parse(Int,ARGS[11]),parse(Bool,ARGS[12]),parse(Bool,ARGS[13]))

# Arguments
# 1: nchains
# 2: gene
# 3: cell
# 4: cond
# 5: G
# 6: maxtime
# 7: infolder
# 8: resultfolder
# 9: datafolder
# 10: inlabel
# 11: label
# 12: nsets (number of rate sets)
# 13: runcycle (bool)
# 14: transient (bool)
