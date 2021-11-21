# Script executed by Biowulf to run code

@everywhere include("/home/carsonc/StochasticGene/src/StochasticGene.jl")

include("/home/carsonc/StochasticGene/runfiles/scriptfunctions.jl")
@time fit_rna(parse(Int,ARGS[1]),ARGS[2],ARGS[3],ARGS[15],ARGS[4],parse(Int,ARGS[5]),parse(Float64,ARGS[6]),ARGS[7],ARGS[8],ARGS[9],ARGS[10],ARGS[11],parse(Int,ARGS[12]),parse(Bool,ARGS[13]),parse(Bool,ARGS[14]))

# fit_rna(nchains::Int,gene::String,cell::String,fittedparam::String,fixedeffects::String,datacond,G::Int,maxtime::Float64,infolder::String,resultfolder::String,datafolder,inlabel,label,nsets,runcycle::Bool=false,transient::Bool=false,samplesteps::Int=100000,warmupsteps=20000,annealsteps=0,temp=1.,tempanneal=100.,root = "/home/carsonc/scrna/")

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
# 15: fittedparams (string), e.g. "1-2-4-7"
