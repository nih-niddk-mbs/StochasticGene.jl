@everywhere include("/home/carsonc/StochasticGene/src/StochasticGene.jl")

include("scriptfunctions.jl")

fit_rna(parse(Int,ARGS[12]),ARGS[1],ARGS[2],parse(Int,ARGS[3]),parse(Float64,ARGS[4]),ARGS[5],ARGS[6],ARGS[7],ARGS[8],ARGS[9],parse(Int,ARGS[10]),parse(Bool,ARGS[11]),parse(Int,ARGS[13]),parse(Int,ARGS[14]),parse(Int,ARGS[15]),parse(Float64,ARGS[16]),parse(Float64,ARGS[17]))
