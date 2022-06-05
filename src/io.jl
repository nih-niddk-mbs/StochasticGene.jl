# io.jl
### Files for saving and reading mh results

abstract type Fields end

struct Result_Fields <: Fields
    name::String
    label::String
    cond::String
    gene::String
    model::String
    nalleles::String
end

struct Summary_Fields <: Fields
    name::String
    label::String
    cond::String
    model::String
end

const raterow_dict = Dict([("ml", 1),("mean",2),("median",3),("last",4)])
const statrow_dict = Dict([("mean",1),("SD", 2),("median",3),("MAD",4)])

"""
  write_dataframes(resultfolder::String,datafolder::String;measure::Symbol=:AIC,assemble::Bool=true)

  collates run results into a csv file

Arguments
- `resultfolder`: name of folder with result files
- `datafolder`: name of folder where data is stored
- `measure`: measure used to assess winner
- `assemble`: if true then assemble results into summary files
"""



function write_dataframes(resultfolder::String,datafolder::String;measure::Symbol=:AIC,assemble::Bool=true)
    write_dataframes_only(resultfolder,datafolder,assemble=assemble)
    write_winners(resultfolder,measure)
end

function write_dataframes_only(resultfolder::String,datafolder::String;assemble::Bool=true)
    dfs = make_dataframes(resultfolder,datafolder,assemble)
    for df in dfs
        for dff in dfs
            for dfff in dff
                csvfile = joinpath(resultfolder,dfff[1])
                CSV.write(csvfile,dfff[2])
            end
        end
    end
    nothing
end

"""
write_winners(resultfolder,measure)

Write best performing model for measure

"""
function write_winners(resultfolder,measure)
    df = best_measure(resultfolder,measure)
    for i in eachindex(df)
        csvfile = joinpath(resultfolder,df[i][1])
        CSV.write(csvfile,df[i][2])
    end
    nothing
end

"""
write_augmented(summaryfile::String,resultfolder,datafolder;fishdata=false)

Augment summary file with G=2 burst size, model predicted moments, and fit measures


"""
write_augmented(summaryfile::String,resultfolder;fishdata=false) = CSV.write(summaryfile,augment_dataframe(read_dataframe(summaryfile),resultfolder))

"""
read_dataframe(csvfile::String)
"""
read_dataframe(csvfile::String) = DataFrame(CSV.File(csvfile))

"""
get_suffix(file::String)

"""
get_suffix(file::String) = chop(file,tail=4), last(file,3)

# does not account for csv files with less than 4 fields
function fields(file::String)
    file,suffix = get_suffix(file)
    v = split(file,"_")
    if suffix == "csv"
        if length(v) == 4
            s = Summary_Fields(v[1],v[2],v[3],v[4])
        else
            println(file)
            throw("Incorrect file name format")
        end
    else
        if length(v) == 6
            s = Result_Fields(v[1],v[2],v[3],v[4],v[5],v[6])
        else
            println(file)
            throw("Incorrect file name format")
        end
    end
    return s
end

function isratefile(folder::String)
    files=readdir(folder)
    any(occursin.(".csv",files) .& occursin.("rates",files))
end

isfish(string::String) = occursin("FISH",string)

function get_genes(file::String)
    r,header = readdlm(file,',',header=true)
    return r[:,1]
end

get_genes(root,cond,datafolder) = get_genes(cond,joinpath(root,datafolder))

function get_genes(cond,datafolder)
    genes = Vector{String}(undef,0)
    files = readdir(datafolder)
    for file in files
        if occursin(cond,file)
            push!(genes,split(file,"_")[1])
        end
    end
    return genes
end
"""
    get_genes(folder,type,label,cond,model)

"""
# not working
function get_genes(folder,type,label,cond,model)
    genes = Array{String,1}(undef,0)
    files = get_files(folder,type,label,cond,model)
    for file in files
        push!(genes,get_gene(file))
    end
    return genes
end

get_files(folder::String,resultname,label,cond,model) = get_files(get_resultfiles(folder),resultname,label,cond,model)

function get_files(files::Vector,resultname,label,cond,model)
    parts = fields.(files)
    files[(getfield.(parts,:name) .== resultname) .& (getfield.(parts,:label) .== label) .& (getfield.(parts,:cond) .== cond) .& (getfield.(parts,:model) .== model)]
end

get_gene(file::String) = fields(file).gene
get_model(file::String) = fields(file).model
get_label(file::String) = fields(file).label
get_cond(file::String) = fields(file).cond
get_nalleles(file::String) = fields(file).nalleles

get_fields(parts::Vector{T}, field::Symbol) where T <: Fields = unique(getfield.(parts,field))

get_models(parts::Vector{T}) where T <: Fields = get_fields(parts,:model)

get_genes(parts::Vector{T}) where T <: Fields = get_fields(parts,:gene)

get_conds(parts::Vector{T}) where T <: Fields = get_fields(parts,:cond)

get_labels(parts::Vector{T}) where T <: Fields = get_fields(parts,:label)

get_names(parts::Vector{T}) where T <: Fields = get_fields(parts,:name)

get_nalleles(parts::Vector{T}) where T <: Fields = get_fields(parts,:nalleles)

get_resultfiles(folder::String) = get_resultfiles(readdir(folder))
get_resultfiles(files::Vector) = files[occursin.(".txt",files) .& occursin.("_",files)]

get_summaryfiles(folder::String) = get_summaryfiles(readdir(folder))
get_summaryfiles(files::Vector) = files[occursin.(".csv",files) .& occursin.("_",files)]
get_summaryfiles(files::Vector,name) = files[occursin.(".csv",files) .& occursin.(name,files)]

get_ratesummaryfiles(files::Vector) = get_summaryfiles(files,"rates")
get_ratesummaryfiles(folder::String) = get_ratesummaryfiles(get_summaryfiles(folder))

get_measuresummaryfiles(files::Vector) = get_summaryfiles(files,"measures")
get_measuresummaryfiles(folder::String) = get_measuresummaryfiles(get_summaryfiles(folder))

get_burstsummaryfiles(files::Vector) = get_summaryfiles(files,"burst")
get_burstsummaryfiles(folder::String) = get_burstsummaryfiles(get_summaryfiles(folder))

"""
write_moments(outfile,genelist,cond,datafolder,fish,root)

"""
function write_moments(outfile,genelist,cond,datafolder,fish,root)
    f = open(outfile,"w")
    writedlm(f,["Gene" "Expression Mean" "Expression Variance"],',')
    for gene in genelist
        h = get_histogram_rna(gene,cond,datafolder,fish,root)
        writedlm(f,[gene mean_histogram(h) var_histogram(h)],',')
    end
    close(f)
end

"""
    write_burst_stats(outfile,infile::String,G::String,cell,folder,cond,root)

"""
function write_burst_stats(outfile,infile::String,G::String,cell,folder,cond,root)
    folder = joinpath(root,folder)
    condarray = split(cond,"-")
    g = parse(Int,G)
    lr = 2*g
    lc = 2*g-1
    freq = Array{Float64,1}(undef,2*length(condarray))
    burst = similar(freq)
    f = open(joinpath(folder,outfile),"w")
    contents,head = readdlm(joinpath(folder,infile),',',header=true)
    label = Array{String,1}(undef,0)
    for c in condarray
        label = vcat(label, "Freq " * c, "sd","Burst Size " * c, "sd")
    end
    writedlm(f,["gene" reshape(label,1,length(label))],',')
    for r in eachrow(contents)
        gene = String(r[1])
        rates = r[2:end]
        rdecay = decay(root,cell,gene)
        cov = read_covparam(joinpath(folder,getfile("param-stats",gene,G,folder,cond)[1]))
        # mu = readmean(joinpath(folder,getfile("param-stats",gene,G,folder,cond)[1]))
        if size(cov,2) < 2
            println(gene)
        end
        for i in eachindex(condarray)
            j = i-1
            freq[2*i-1], freq[2*i] = frequency(rates[1+lr*(i-1)],sqrt(cov[1+lc*j,1+lc*j]),rdecay)
            burst[2*i-1], burst[2*i] = burstsize(rates[3+lr*j],rates[2+lr*j],cov[3+lc*j,3+lc*j],cov[2+lc*j,2+lc*j],cov[2+lc*j,3+lc*j])
        end
        writedlm(f,[gene freq[1] freq[2] burst[1] burst[2] freq[3] freq[4] burst[3] burst[4]],',')
        flush(f)
    end
    close(f)
end
"""
    write_histograms(resultfolder,ratefile,cell,datacond,G::Int,datafolder::String,fish,root,outfolder = "histograms")

"""
function write_histograms(resultfolder,ratefile,cell,datacond,G::Int,datafolder::String,fish,root,outfolder = "histograms")
    ratefile = joinpath(resultfolder,ratefile)
    rates,head = readdlm(ratefile,',',header=true)
    outfolder = joinpath(resultfolder,outfolder)
    if ~isdir(outfolder)
        mkpath(outfolder)
    end
    cond = string.(split(datacond,"-"))
    for r in eachrow(rates)
        h = histograms(r,cell,cond,G,datafolder,fish,root)
        for i in eachindex(cond)
            f = open(joinpath(outfolder,string(r[1]) * cond[i] * ".txt"),"w")
            writedlm(f,h[i])
            close(f)
        end
    end
end


"""
    assemble_all(folder)

"""
function assemble_all(folder::String)
    files = get_resultfiles(folder)
    parts = fields.(files)
    labels = get_labels(parts)
    conds = get_conds(parts)
    models = get_models(parts)
    assemble_all(folder,files,labels,conds,models)
end

function assemble_all(folder::String,files::Vector,labels::Vector,conds::Vector,models::Vector)
    for l in labels, c in conds, g in models
        assemble_all(folder,files,l,c,g,isfish(l))
    end
end

function assemble_all(folder::String,files::Vector,label::String,cond::String,model::String,fish::Bool)
    assemble_rates(folder,files,label,cond,model,fish)
    assemble_measures(folder,files,label,cond,model)
    assemble_stat(folder,files,label,cond,model)
    if model == "2"
        assemble_burst_model2(folder,files,label,cond)
    end
end

function assemble_files(folder::String,files::Vector,outfile::String,header,readfunction)
    f = open(outfile,"w")
    writedlm(f,header,',')
    for file in files
        gene = get_gene(file)
        r = readfunction(joinpath(folder, file))
        writedlm(f,[gene r],',')
    end
    close(f)
end

function assemble_rates(folder::String,files::Vector,label::String,cond::String,model::String,fish)
    outfile = joinpath(folder,"rates_" * label * "_" * cond * "_" * model * ".csv")
    assemble_files(folder,get_files(files,"rates",label,cond,model),outfile,ratelabels(model,length(split(cond,"-"))),readml)
end

function assemble_measures(folder::String,files,label::String,cond::String,model::String)
    outfile = joinpath(folder,"measures_" * label * "_" * cond * "_" * model *  ".csv")
    header = ["Gene" "Nalleles" "Deviance" "LogMaxLikelihood" "WAIC" "SD" "AIC" "Acceptance" "Temperature" "Rhat"]
    # assemble_files(folder,get_files(files,"measures",label,cond,model),outfile,header,readmeasures)
    files = get_files(files,"measures",label,cond,model)
    f = open(outfile,"w")
    writedlm(f,header,',')
    for file in files
        gene = get_gene(file)
        nalleles = get_nalleles(file)
        r = readmeasures(joinpath(folder, file))
        writedlm(f,[gene nalleles r],',')
    end
    close(f)
end

function assemble_mad(folder::String,files,label::String,cond::String,model::String)
    outfile = joinpath(folder,"mad_" * label * "_" * cond * "_" * model *  ".csv")
    assemble_files(folder,get_files(files,"param-stats",label,cond,model),outfile,ratelabels(model,length(split(cond,"-")))[:,1:end-1],readmad)
end

function assemble_sd(folder::String,files,label::String,cond::String,model::String)
    outfile = joinpath(folder,"sd_" * label * "_" * cond * "_" * model *  ".csv")
    assemble_files(folder,get_files(files,"param-stats",label,cond,model),outfile,ratelabels(model,length(split(cond,"-")),"SD")[:,1:end-1],readsd)
end

function assemble_stat(folder::String,files,label::String,cond::String,model::String)
    outfile = joinpath(folder,"stats_" * label * "_" * cond * "_" * model *  ".csv")
    assemble_files(folder,get_files(files,"param-stats",label,cond,model),outfile,statlabels(model,length(split(cond,"-"))),readstats)
end

function assemble_burst_model2(folder::String,files,label::String,cond::String,model::String = "2")
    outfile = joinpath(folder,"burst_" * label * "_" * cond * "_" * model *  ".csv")
    assemble_files(folder,get_files(files,"param-stats",label,cond,model),outfile,["Gene" "BurstSize" "BurstSD" "BurstVar" "E" "O" "VEE" "VOO" "VEO"],read_burst_model2)
end

function Gratelabels(model,nsets)
    G = parse(Int,model)
    n = G-1
    Grates = Array{String,2}(undef,1,2*n)
    for i = 0:n-1
        Grates[1,2*i+1] = "Rate$i$(i+1)"
        Grates[1,2*i+2] = "Rate$(i+1)$i"
    end
    return Grates
end

function ratelabels(model,nsets)
    Grates = Gratelabels(model,nsets)
    rates = [Grates "Eject" "Decay"]
    for i = 2:nsets
        rates = [rates rates]
    end
    return ["Gene" rates]
end

function statlabels(model,nsets)
    label = ["Mean","SD","Median","MAD"]
    Grates = Gratelabels(model,nsets)
    Grates = [Grates "Eject"]
    for i = 2:nsets
        Grates = [Grates Grates]
    end
    rates = Matrix{String}(undef,1,0)
    for i in 1:4
        rates = [rates Grates .* label[i]]
    end
    return ["Gene" rates]
end

# function ratelabels(model,nsets,sd::Bool=false)
#     G = parse(Int,model)
#     n = G-1
#     Grates = Array{String,2}(undef,1,2*n)
#     for i = 0:n-1
#         Grates[1,2*i+1] = "Rate$i$(i+1)"
#         Grates[1,2*i+2] = "Rate$(i+1)$i"
#     end
#     rates = [Grates "Eject" "Decay"]
#     for i = 2:nsets
#         rates = [rates Grates "Eject" "Decay"]
#     end
#     if typeof(model) <: AbstractGMlossmodel
#         rates = [rates "Yield"]
#         lenr = 2*G*nsets + 1
#     else
#         lenr = 2*G*nsets
#     end
#     if sd
#         rates = reshape([rates; fill("SD",1,lenr)],1,2*(lenr))
#     end
#     return ["Gene" rates]
# end

function get_all_rates(file::String,header::Bool)
    r = readdlm(file,',',header=header)
    if header
        r = r[1]
    end
    return r
end



filename(data,model::AbstractGRMmodel) = filename(data.name,data.gene,model.G,model.R,model.nalleles)
filename(data,model::AbstractGMmodel) = filename(data.name,data.gene,model.G,model.nalleles)
filename(label::String,gene::String,G::Int,R::Int,nalleles::Int) = filename(label,gene,"$G"*"$R","$(nalleles)")
filename(label::String,gene,G::Int,nalleles::Int) = filename(label,gene,"$G","$(nalleles)")
filename(label::String,gene::String,model::String,nalleles::String) = "_" * label  * "_" * gene *  "_" * model * "_" * nalleles * txtstr


"""
write_results(file::String,x)
"""
function writeall(path::String,fit,stats,measures,data,temp,model::StochasticGRmodel)
    if ~isdir(path)
        mkpath(path)
    end
    name = filename(data,model)
    write_rates(joinpath(path,"rates" * name ),fit,stats,model)
    write_measures(joinpath(path,"measures" * name),fit,measures,deviance(fit,data,model),temp)
    write_param_stats(joinpath(path,"param-stats" * name),stats)

end

"""
write_rates(file::String,fit)

Write rate parameters, rows in order are
maximum likelihood
mean
median
last accepted
"""
function write_rates(file::String,fit::Fit,stats,model)
    f = open(file,"w")
    writedlm(f,[get_rates(fit.parml,model)],',')
    writedlm(f,[get_rates(stats.meanparam,model)],',')
    writedlm(f,[get_rates(stats.medparam,model)],',')
    writedlm(f,[get_rates(fit.param[:,end],model)],',')
    close(f)
end
"""
write_measures(file,fit,waic,dev)
"""
function write_measures(file::String,fit::Fit,measures::Measures,dev,temp)
    f = open(file,"w")
    writedlm(f,[fit.llml mean(fit.ll) std(fit.ll) quantile(fit.ll,[.025;.5;.975])' measures.waic[1] measures.waic[2] aic(fit)],',')
    writedlm(f,dev,',')
    writedlm(f,[fit.accept fit.total],',')
    writedlm(f,temp,',')
    writedlm(f,measures.rhat',',')
    writedlm(f,maximum(measures.rhat),',')
    close(f)
end
"""
write_param_stats(stats,waic,data,model)

"""
function write_param_stats(file,stats::Stats)
    f = open(file,"w")
    writedlm(f,stats.meanparam',',')
    writedlm(f,stats.stdparam',',')
    writedlm(f,stats.medparam',',')
    writedlm(f,stats.madparam',',')
    writedlm(f,stats.qparam,',')
    writedlm(f,stats.corparam,',')
    writedlm(f,stats.covparam,',')
    writedlm(f,stats.covlogparam,',')
    close(f)
end

"""
readrates(file::String,row::Int)

row
1       maximum likelihood
2       mean
3       median
4       last value of previous run
"""

readrates(file::String) = readrates(file,3)
readrates(file::String,row::Int) = readrow(file,row)

function readrow(file::String,row)
    if isfile(file) && ~isempty(read(file))
        contents = readdlm(file,',')
        if row <= size(contents,1)
            m = contents[row,:]
            return m[.~isempty.(m)]
        else
            return [nothing]
        end
    else
        println(file)
    end
end

function readrow_flip(file,row)
    m = readrow(file,row)
    reshape(m,1,length(m))
end

function readmeasures(file::String)
    d = readdeviance(file)
    w = readwaic(file)
    a = readaccept(file)
    t = readtemp(file)
    r = readrhat(file)
    [d[1] w[1] w[7] w[8] w[9] a t[1] r[1]]
end

readdeviance(file::String) = readrow(file,2)

readwaic(file::String) = readrow(file,1)

function readaccept(file::String)
    a = readrow(file,3)
    a[1]/a[2]
end

readtemp(file::String) = readrow(file,4)

readrhat(file::String) = readrow(file,6)

function readml(ratefile::String)
    m = readrow(ratefile,1)
    reshape(m,1,length(m))
end

function readmean(ratefile::String)
    m = readrow(ratefile,2)
    reshape(m,1,length(m))
end

function readsd(ratefile::String)
    m = readrow(ratefile,2)
    reshape(m,1,length(m))
end

function readstats(statfile::String)
    mean = readrow_flip(statfile,1)
    sd = readrow_flip(statfile,2)
    median = readrow_flip(statfile,3)
    mad = readrow_flip(statfile,4)
    [mean sd median mad]
end

function readmedian(statfile::String)
    m = readrow(statfile,3)
    reshape(m,1,length(m))
end

function readmad(statfile::String)
    m = readrow(file,4)
    reshape(m,1,length(m))
end

function read_corparam(file::String)
    c = readdlm(file,',')
    n = length(c[1,:])
    # c[5+n:4+2*n,1:n]
    c[8:7+n,1:n]
end

function read_covparam(file::String)
    c = readdlm(file,',')
    read_covparam(c)
end

function read_covparam(c::Matrix)
    n = length(c[1,:])
    # c[5+2*n:4+3*n,1:n]
    c[8+n:7+2*n,1:n]
end

function read_covlogparam(file::String)
    c = readdlm(file,',')
    n = length(c[1,:])
    c[8+2*n:7+3*n,1:n]
end

read_crosscov(statfile::String) = read_crosscov(read_covparam(statfile))

function read_crosscov(C::Matrix)
    c = Float64[]
    N = size(C,1)
    for i in 1:N
        for j in i+1:N
            push!(c,C[i,j])
        end
    end
    c
end


function read_burst_model2(file::String)
    c = readdlm(file,',')
    b = c[1,end]/c[1,end-1]
    cov = read_covparam(c)
    v = var_ratio(c[1,end],c[1,end-1],cov[end,end],cov[end-1,end-1],cov[end-1,end])
    [b sqrt(abs(v)) v c[1,end] c[1,end-1] cov[end,end] cov[end-1,end-1] cov[end-1,end]]
end

function write_residency_G(fileout::String,filein::String,G,header)
    rates = get_all_rates(filein,header)
    n = G-1
    f = open(fileout,"w")
    writedlm(f,["gene" collect(0:n)'],',')
    for r in eachrow(rates)
        writedlm(f,[r[1] residenceprob_G(r[2:2*n+1],n)],',')
    end
    close(f)
end


# Functions for saving and loading data and models

"""
write_log(file,datafile,data,model)
write all information necessary for rerunning
"""
function save_data(file::String,data::TransientRNAData)
    f = open(file,"w")
    writedlm(f, [typeof(data)])
    writedlm(f,[data.name])
    writedlm(f,[data.gene])
    writedlm(f,[data.nRNA])
    writedlm(f,[data.time])
    writedlm(f,[data.histRNA])
    close(f)
end

function load_data(file::String,model::AbstractGMmodel)


end

function save_model(file::String,model::AbstractGMmodel)
    f = open(file,"w")
    write(f, model.G)
    writedlm(f,model.nalleles)
    writedlm(f,model.ratepriors)
    writedlm(f,model.proposal)
    writedlm(f,model.fittedparam)
    writedlm(f,model.method)
    close(f)

end

function load_model(file::String,model::AbstractGRMmodel)

end

"""
    makestring(v::Vector)

    turn a vector of strings into a single string
"""
function makestring(v)
    s = ""
    for i in v
        s *= i
    end
    return s
end
