# Functions for creating batch files to run on NIH Biowulf super computer

using DelimitedFiles

transientgenelist = ["MYC", "DUSP5", "TRIB1", "PMAIP1", "SERPINE1", "SOX9", "ERRFI1", "NR4A2", "JUN", "GEM"]


function makeswarm(G::Int,infolder,swarmfile::String,inlabel,label,nsets,datafolder::String,thresholdlow,thresholdhigh,conds::Vector = ["DMSO","AUXIN"],result::String= "2021-03-11",batchsize=1000,maxtime = 3600. * 8,nchains::Int = 8,runcycle::Bool=true,transient::Bool=false,juliafile::String="/home/carsonc/StochasticGene/runfiles/fitscript.jl",nprocs::Int=8,root="/home/carsonc/scrna/")
    genes = checkgenes(conds[1],datafolder,thresholdlow,thresholdhigh,root)
    makeswarm(genes,G,infolder,swarmfile,inlabel,label,nsets,datafolder,conds,result,batchsize,maxtime,nchains,runcycle,transient,juliafile,nprocs,root)
end

function makeswarm(genes,G::Int,infolder,swarmfile::String,inlabel,label,nsets,datafolder::String,conds::Vector,result::String,batchsize,maxtime,nchains::Int,runcycle::Bool=true,transient::Bool=false,juliafile::String="fitscript.jl",nprocs::Int=8,root="/home/carsonc/scrna/")

    resultfolder = joinpath("Results",result)
    infolder = joinpath("Results",infolder)
    ngenes = length(genes)
    println(ngenes)
    println(runcycle)

    if ngenes > batchsize
        nbatches = div(ngenes,batchsize)
        batches = Vector{Vector{String}}(undef,nbatches+1)
        println(batchsize," ",nbatches)
        for i in 1:nbatches
            batches[i] = genes[batchsize*(i-1)+1:batchsize*(i)]
        end
        batches[end] = genes[batchsize*nbatches+1:end]
        for cond in conds
            for batch in eachindex(batches)
                sfile = swarmfile * "_" * cond * "$batch" * ".swarm"
                writegenes(sfile,batches[batch],nprocs,juliafile,cond,G,maxtime,infolder,resultfolder,datafolder,inlabel,label,nsets,runcycle,transient,nchains)
            end
        end
    else
        for cond in conds
            sfile = swarmfile * "_" * cond * ".swarm"
            f = open(sfile,"w")
            writegenes(sfile,genes,nprocs,juliafile,cond,G,maxtime,infolder,resultfolder,datafolder,inlabel,label,nsets,runcycle,transient,nchains)
        end
    end
end

# function makeswarm(G::Int,infolder,swarmfile::String,dataname::String,datafolder::String,result::String,batchsize::Int=1000,maxtime = 3600. * 8,nchains = 8,runcycle::Bool=true,transient::Bool,juliafile::String="fitscript.jl",nprocs::Int=8,root="/home/carsonc/scrna/")
#     genes = checkgenes(datafolder,root)
#     makeswarm(genes,G,infolder,swarmfile,dataname,datafolder,result,batchsize,maxtime,nchains,runcycle,transient,juliafile,nprocs::Int,root)
# end
#
# function makeswarm(genes,G::Int,infolder::String,swarmfile::String,dataname::String,datafolder::String,result::String,batchsize,maxtime,nchains::Int,runcycle::Bool,transient::Bool,juliafile::String="fitscript.jl",nprocs::Int=8,root::String="/home/carsonc/scrna/")
#
#     # genes = checkgenes(datafolder,root)
#     resultfolder = joinpath("Results",result)
#     infolder = joinpath("Results",infolder)
#     label = "scRNA_" * dataname * "_ss"
#     inlabel = "scRNA_" * dataname * "_ss"
#     nsets = 1
#
#     ngenes = length(genes)
#     println(ngenes)
#     println(runcycle)
#
#     if ngenes > batchsize
#         batches = getbatches(genes,ngenes,batchsize)
#         for batch in eachindex(batches)
#             sfile = swarmfile * "_" * "$batch" * ".swarm"
#             writegenes(sfile,batches[batch],nprocs,juliafile,"null",G,maxtime,infolder,resultfolder,datafolder,inlabel,label,nsets,runcycle,transient,nchains)
#         end
#     else
#         sfile = swarmfile * ".swarm"
#         writegenes(sfile,genes,nprocs,juliafile,"null",G,maxtime,infolder,resultfolder,datafolder,inlabel,label,nsets,runcycle,transient,nchains)
#     end
# end

function getbatches(genes,ngenes,batchsize)
    nbatches = div(ngenes,batchsize)
    batches = Vector{Vector{String}}(undef,nbatches+1)
    println(batchsize," ",nbatches)
    for i in 1:nbatches
        batches[i] = genes[batchsize*(i-1)+1:batchsize*(i)]
    end
    batches[end] = genes[batchsize*nbatches+1:end]
    return batches
end

function writegenes(sfile,genes,nprocs,juliafile,cond,G,maxtime,infolder,resultfolder,datafolder,inlabel,label,nsets,runcycle,transient,nchains)
    f = open(sfile,"w")
println(runcycle)
    for gene in genes
        writedlm(f,["julia -p" nprocs juliafile gene cond G maxtime infolder resultfolder datafolder inlabel label nsets runcycle transient nchains])
    end
    close(f)
end

function checkgenes(cond,datafolder,thresholdlow,thresholdhigh,root = "/Users/carsonc/Box/scrna")
    datapath = joinpath(root,datafolder)
    # genes = readdlm(joinpath(root,"data/HCT_HEK_Half_life_top_below_1h.csv"),',')[2:end,end]
    genedecay = readdlm(joinpath(root,"data/HCT116_all_cells_histograms_and_half_lives_March_2021/Starved_Rep7_half_lives.csv"),',')[2:end,2]
    genenames = readdlm(joinpath(root,"data/HCT116_all_cells_histograms_and_half_lives_March_2021/Starved_Rep7_half_lives.csv"),',')[2:end,1]
    genes1 = Vector{String}(undef,0)
    for i in eachindex(genedecay)
        if typeof(genedecay[i]) <: Number
            if thresholdlow <= genedecay[i] < thresholdhigh
                push!(genes1,genenames[i])
            end
        end
    end
    println(length(genes1))
    genes2 = readdlm(joinpath(root,"data/HCT116_alleles_number.txt"))[2:end,1]
    println(length(genes2))
    genes = intersect(genes1,genes2)
    complement = setdiff(genes1,genes2)
    println(length(genes))
    println(length(complement))
    list = Vector{String}(undef,0)
    for gene in genes
        if cond == ""
            datafile = joinpath(datapath,gene * ".txt")
        else
            datafile = joinpath(datapath,gene * "_" * cond * ".txt")
        end
        if isfile(datafile)
            push!(list,gene)
        end
    end
    println(length(list))
    return list
end


function checkgenes(cond,datafolder,root = "/Users/carsonc/Box/scrna")
    datapath = joinpath(root,datafolder)
    # genes = readdlm(joinpath(root,"data/HCT_HEK_Half_life_top_below_1h.csv"),',')[2:end,end]
    genedecay = readdlm(joinpath(root,"data/HCT116_all_cells_histograms_and_half_lives_March_2021/Starved_Rep7_half_lives.csv"),',')
    # genes1 = readdlm(joinpath(root,"data/HCT116_all_cells_histograms_and_half_lives_March_2021/Starved_Rep7_half_lives.csv"),',')[2:end,1]
    genes1 = Vector{String}(undef,0)
    for i in eachindex(genedecay[2:end,2])
        if typeof(genedecay[i,2]) <: Number
            push!(genes1,genedecay[i,1])
        end
    end
    genes2 = readdlm(joinpath(root,"data/HCT116_alleles_number.txt"))[2:end,1]
    println(length(genes1))
    println(length(genes2))
    genes = intersect(genes1,genes2)
    list = Vector{String}(undef,0)
    for gene in genes
        datafile = joinpath(datapath,gene * "_" * cond * ".txt")
        if isfile(datafile)
            push!(list,gene)
        end
    end
    println(length(list))
    return list
end


function checkgenes(datafolder,root)
    datapath = joinpath(root,datafolder)
    genedecay = readdlm(joinpath(root,"data/HCT116_all_cells_histograms_and_half_lives_March_2021/Starved_Rep7_half_lives.csv"),',')
    genes1 = Vector{String}(undef,0)
    for i in eachindex(genedecay[2:end,2])
        if typeof(genedecay[i,2]) <: Number
            push!(genes1,genedecay[i,1])
        end
    end
    genes2 = readdlm(joinpath(root,"data/HCT116_alleles_number.txt"))[2:end,1]
    println(length(genes1))
    println(length(genes2))
    genes = intersect(genes1,genes2)
    list = Vector{String}(undef,0)
    for gene in genes
        datafile = joinpath(datapath,gene * ".txt")
        if isfile(datafile)
            push!(list,gene)
        end
    end
    println(length(list))
    return list
end

fix(folder) = writeruns(fixruns(findjobs(folder)))


function findjobs(folder)
    files = readdir(folder)
    files = files[occursin.("swarm_",files)]
    for (i,file) in enumerate(files)
        files[i] = split(file,"_")[2]
    end
    unique(files)
end

function fixruns(jobs)
    runlist = Vector{String}(undef,0)
    for job in jobs
        if occursin("FAILED",read(`jobhist $job`,String))
            swarmfile = findswarm(job)
            list = readdlm(swarmfile,',')
            runs =  chomp(read(pipeline(`jobhist $job`, `grep FAILED`),String))
            runs = split(runs,'\n')
            println(job)
            for run in runs
                linenumber = parse(Int,split(split(run," ")[1],"_")[2])
                a = String(list[linenumber+1])
                println(a)
                push!(runlist,a)
            end
        end
    end
    return runlist
end

function writeruns(runs)

    f = open("fitfix.swarm","w")
    for run in runs
        writedlm(f,[run],quotes=false)
    end
    close(f)

end

function findswarm(job)
    sc = "Swarm Command"
    line = read(pipeline(`jobhist $job`, `grep $sc`),String)
    list = split(line," ")
    list[occursin.(".swarm",list)][1]
end

function collate(folders,type="rates_scRNA_T120")
    genes = Vector{Vector}(undef,0)
    for folder in folders
        files = readdir(folder)
        for file in files
            if occursin(type,file)
                contents=readdlm(file,',')
                for row in eachrow(contents)
                    row[1] .== genes
                end
            end
        end
    end
end

large_deviance(measurefile,threshold) = filter_gene(measurefile,"Deviance",threshold)

function filter_gene(measurefile,measure,threshold)
    genes = Vector{String}(undef,0)
    measures,header = readdlm(measurefile,',',header=true)
    println(length(measures[:,1]))
    col = findfirst(header[1,:] .== measure)
    for d in eachrow(measures)
        if d[col] > threshold || isnan(d[col])
            push!(genes,d[1])
        end
    end
    println(length(genes))
    return genes
end

function filter_gene_nan(measurefile,measure)
    genes = Vector{String}(undef,0)
    measures,header = readdlm(measurefile,',',header=true)
    println(length(measures[:,1]))
    col = findfirst(header[1,:] .== measure)
    for d in eachrow(measures)
        if isnan(d[col])
            push!(genes,d[1])
        end
    end
    println(length(genes))
    return genes
end


function replace_yield(G,folder1,folder2,cond1,cond2,outfolder)
    if typeof(G) <: Number
        G = string(G)
    end
    if ~isdir(outfolder)
        mkpath(outfolder)
    end
    files1 = getratefile(folder1,G,cond1)
    files2 = getratefile(folder2,G,cond2)
    for file1 in files1
        gene = StochasticGene.getgene(file1)
        file2 = getratefile(files2,gene)
        outfile = joinpath(outfolder,file2)
        r1 = StochasticGene.readrates(joinpath(folder1,file1))
        r2 = StochasticGene.readrates(joinpath(folder2,file2))
        r2[end] = r1[end]
        f = open(outfile,"w")
        writedlm(f,[r2],',')
        close(f)
    end

end

"""
assemble_r(G,folder1,folder2,cond1,cond2,outfolder)

Combine rates from two separate fits into a single rate vector

"""

function assemble_r(G,folder1,folder2,cond1,cond2,outfolder)
    if typeof(G) <: Number
        G = string(G)
    end
    if ~isdir(outfolder)
        mkpath(outfolder)
    end
    files1 = getratefile(folder1,G,cond1)
    files2 = getratefile(folder2,G,cond2)
    for file1 in files1
        gene = StochasticGene.getgene(file1)
        file2 = getratefile(files2,gene)
        name = replace(file1, cond1 => "JOINT")
        outfile = joinpath(outfolder,name)
        assemble_r(joinpath(folder1,file1),joinpath(folder2,file2),outfile)
    end
end


function  assemble_r(ratefile1,ratefile2,outfile)
    r1 = StochasticGene.readrates(ratefile1)
    r2 = StochasticGene.readrates(ratefile2)
    r1[end] = clamp(r1[end],eps(Float64),1-eps(Float64))
    r = vcat(r1[1:end-1],r2[1:end-1],r1[end])
    f = open(outfile,"w")
    writedlm(f,[r],',')
    close(f)
end


function getratefile(files,gene)
    files = files[occursin.("_"*gene*"_",files)][1]
end

function getratefile(folder,G,cond)
    files = readdir(folder)
    files = files[occursin.("rates_",files)]
    files = files[occursin.("_"*cond*"_",files)]
    files[occursin.("_"*G*"_",files)]
end

function getratefile(gene,G,folder,cond)
    files = readdir(folder)
    files = files[occursin.("rate",files)]
    files = files[occursin.("_"*gene*"_",files)]
    files = files[occursin.("_"*G*"_",files)]
    files[occursin.(cond,files)]
end



function assemble_r(gene,G,folder1,folder2,cond1,cond2,outfolder)
    file1 = getratefile(gene,G,folder1,cond1)[1]
    file2 = getratefile(gene,G,folder2,cond2)[1]
    name = replace(file1, cond1 => "JOINT")
    println(name)
    outfile = joinpath(outfolder,name)
    println(outfile)
    assemble_r(joinpath(folder1,file1),joinpath(folder2,file2),outfile)
end
