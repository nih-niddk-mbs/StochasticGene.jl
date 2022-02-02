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

"""
  write_dataframes(resultfolder::String,datafolder::String)

  collates run results into a csv file

Arguments
- `resultfolder`: name of folder with result files
- `datafolder`: name of folder where data is stored
"""

function write_dataframes(resultfolder::String,datafolder::String,measure=:AIC)
    dfs = make_dataframes(resultfolder,datafolder)
    for df in dfs
        csvfile = joinpath(resultfolder,df[1])
        CSV.write(csvfile,df[2])
    end
    df = best_measure(resultfolder,measure)
    for i in eachindex(df)
        csvfile = joinpath(resultfolder,df[i][1])
        CSV.write(csvfile,df[i][2])
    end
    nothing
end


function fields(file::String)
    file,suffix = split(file,".")
    v = split(file,"_")
    if suffix == "csv"
        if length(v) == 4
            s = Summary_Fields(v[1],v[2],v[3],v[4])
        else
            throw("Incorrect file name format")
        end
    else
        if length(v) == 6
            s = Result_Fields(v[1],v[2],v[3],v[4],v[5],v[6])
        else
            throw("Incorrect file name format")
        end
    end
    return s
end

isfish(string::String) = occursin("FISH",string)

# function extractparts(file::String)
#     file  = split(file,".")[1]
#     split(file,"_")
# end


function get_genes(folder,type,label,cond,model)
    genes = Array{String,1}(undef,0)
    files = get_files(folder,type,label,cond,model)
    for file in files
        push!(genes,get_gene(file))

    end
    return genes
end

# get_attributes(parts::Vector,f) = unique(map(x->f(x),parts))

get_gene(file::String) = fields(file).gene
get_model(file::String) = fields(file).model

get_fields(parts::Vector{T}, field::Symbol) where T <: Fields = unique(getfield.(parts,field))

get_models(parts::Vector{T}) where T <: Fields = get_fields(parts,:model)

get_genes(parts::Vector{T}) where T <: Fields = get_fields(parts,:gene)

get_conds(parts::Vector{T}) where T <: Fields = get_fields(parts,:cond)

get_labels(parts::Vector{T}) where T <: Fields = get_fields(parts,:label)

get_names(parts::Vector{T}) where T <: Fields = get_fields(parts,:name)

get_files(folder::String,resultname,label,cond,model) = get_files(readdir(folder),resultname,label,cond,model)
function get_files(files::Vector,resultname,label,cond,model)
    parts = fields.(files)
    files[(getfield.(parts,:name) .== resultname) .& (getfield.(parts,:label) .== label) .& (getfield.(parts,:cond) .== cond) .& (getfield.(parts,:model) .== model)]
end

get_resultfiles(folder::String) = get_resultfiles(readdir(folder))
get_resultfiles(files::Vector) = files[occursin.(".txt",files) .& occursin.("_",files)]

get_summaryfiles(folder::String) = get_summaryfiles(readdir(folder))
get_summaryfiles(files::Vector) = files[occursin.(".csv",files) .& occursin.("_",files)]
get_summaryfiles(files::Vector,name) = files[occursin.(".csv",files) .& occursin.(name,files)]

get_ratesummaryfiles(files::Vector) = get_summaryfiles(files,"rates")
get_ratesummaryfiles(folder::String) = get_ratesummaryfiles(get_summaryfiles(folder))

get_measuresummaryfiles(files::Vector) = get_summaryfiles(files,"measures")
get_measuresummaryfiles(folder::String) = get_measuresummaryfiles(get_summaryfiles(folder))


# function make_dataframe(labels::Vector,models::Vector,conds::String,resultfolder::String,datafolder::String,fish::Bool,G::Int = 2)
#     # folder = joinpath(root,joinpath("results",folder))
#     cond = [conds]
#     dataconds = string.(split(conds,"-"))
#     datafolder = repeat([datafolder],length(dataconds))
#     if isempty(label)
#         assemble_all(resultfolder,cond,models,fish)
#     else
#         assemble_all(resultfolder,label,cond,models,fish)
#     end
#     df = make_dataframe(folder,models::Vector,fish,G)
#     add_mean(df,dataconds,datafolder,fish,root)
# end

# function make_dataframe_2cond(resultfolder,models::Vector,fish::Bool,G::Int = 2)
#     files = readdir(resultfolder)
#     mfiles = Vector{String}(undef,0)
#     rfile = ""
#     for file in files
#         if occursin("measures",file)  && occursin("csv",file) && ~occursin("all",file)
#             push!(mfiles,joinpath(resultfolder,file))
#         elseif occursin("rates",file) && occursin("$(G).csv",file)
#             rfile = joinpath(folder,file)
#         end
#     end
#     r,head = readdlm(rfile,',',header=true)
#     winner = best_model(models,mfiles)
#     cond = [zeros(length(r[:,1]));ones(length(r[:,1]))];
#     namelist = Symbol.(string.(head))
#     rs = [vcat(r[:,[1;2:2*G+1]], r[:,[1;2*G+2:4*G+1]])  cond]
#     df = DataFrame()
#     df[!,namelist[1]] = string.(rs[:,1])
#     for i in 2:2*G+1
#         df[!,namelist[i]] = float.(rs[:,i])
#     end
#     if ~fish
#         rs = [vcat(r[:,[1,2,3,4,5,10]], r[:,[1,6,7,8,9,10]])  cond winner]
#         vcat(r[:,4*G + 2], r[:,4*G + 2])
#         df = [df DataFrame(yield = float.(vcat(r[:,4*G + 2], r[:,4*G + 2])))]
#     end
#     df = [df DataFrame(cond = Int.(rs[:,end]))]
#     df = innerjoin(df, winner, on = :Gene)
# end

function make_dataframes(resultfolder::String,datafolder::String)
    assemble_all(resultfolder)
    files = get_ratesummaryfiles(resultfolder)
    parts = fields.(files)
    models = get_models(parts)
    labels = get_labels(parts)
    df = Vector{Tuple{String,DataFrame}}(undef,0)
    for label in labels
        lfiles = files[occursin.(label,files)]
        dfl = Vector{DataFrame}(undef,0)
        for model in models
            mfiles = files[occursin.(model,lfiles)]
            dfm = Vector{DataFrame}(undef,length(mfiles))
            for i in eachindex(mfiles)
                dfm[i] = make_dataframe(joinpath(resultfolder,mfiles[i]),datafolder,isfish(mfiles[i]))
            end
            push!(df,("Summary_$(label)_$(model).csv",stack_dataframe(dfm)))
        end
    end
    return df
end


function make_dataframe(ratefile::String,datafolder::String,fish::Bool)
    df = read_dataframe(ratefile)
    filename = splitpath(ratefile)[end]
    parts = fields(filename)
    G = parse(Int,parts.model)
    insertcols!(df, :Model => fill(G,size(df,1)))
    df  = stack_dataframe(df,G,parts.cond,fish)
    add_mean(df,datafolder,fish)
end

function add_measures(df,measurefile::String)
    dm = read_dataframe(measurefile)
    hcat(df,dm[:,[:Deviance,:WAIC,:AIC]])
end

function add_mean(df::DataFrame,datafolder,fish::Bool)
    root = string(split(abspath(datafolder),"data")[1])
    m = Vector{Float64}(undef,length(df.Gene))
    i = 1
    for gene in df.Gene
        m[i] = mean_histogram(get_histogram_rna(string(gene),df[i,:Condition],datafolder,fish,root))
        i += 1
    end
    insertcols!(df, :Expression => m)
end

stack_dataframe(df,G,cond,fish) = stack_dataframe(separate_dataframe(df,G,cond,fish))

function stack_dataframe(df2::Vector{DataFrame})
    df = df2[1]
    for i in 2:length(df2)
        df = vcat(df,df2[i])
    end
    return df
end

function separate_dataframe(df,G,cond,fish)
    conds = split(cond,"-")
    nsets = length(conds)
    df2 = Vector{DataFrame}(undef,nsets)
    for i in 1:nsets
        df2[i] = df[:,[1; 2*G*(i-1) + 2 : 2*G*i + 1;end-Int(~fish): end]]
        rename!(x->split(x,"_")[1],df2[i])
        insertcols!(df2[i], :Condition => fill(string(conds[i]),size(df,1)))
    end
    return df2
end

read_dataframe(csvfile::String) = DataFrame(CSV.File(csvfile))

best_AIC(folder::String) = best_measure(folder,:AIC)

best_WAIC(folder::String) = best_measure(folder,:WAIC)


function best_measure(folder::String,measure::Symbol)
    files = get_measuresummaryfiles(folder)
    parts = fields.(files)
    labels = get_labels(parts)
    df = Vector{Tuple{String,DataFrame}}(undef,0)
    for label in labels
        lfiles = files[occursin.(label,files)]
        dm = Vector{DataFrame}(undef,length(lfiles))
        for i in eachindex(files)
            dm[i] = read_dataframe(joinpath(folder,files[i]))
            insertcols!(dm[i], :Model => fill(get_model(files[i]),size(dm[i],1)))
        end
        push!(df,("Winners_$(label).csv",best_measure(dm,measure)))
    end
    return df
end

function best_measure(dfs::Vector,measure::Symbol)
    l = dfs[1][:,measure]
    model = dfs[1][1, :Model]
    for i in 2:length(dfs)
        l = [l dfs[i][:,measure]]
        model = [model dfs[i][1,:Model]]
    end
    ind = argmin.(eachrow(l))
    DataFrame(Gene = dfs[1].Gene, Winner = model[ind])
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
    # assemble_stats("mean",folder,l,c,g,append,header)
end

function assemble_files(folder::String,files::Vector,outfile::String,header,readfunction)
    f = open(outfile,"w")
    writedlm(f,header,',')
    for file in files
        gene = get_gene(file)
        r = readfunction(joinpath(folder, file))
        writedlm(f,[gene r'],',')
    end
    close(f)
end

function assemble_files(folder::String,files::Vector,header,readfunction)
    namelist = Symbol.(header)
    df = DataFrame()
       df[!,namelist[1]] = String[]
       for i in 2:length(header)
           df[!,namelist[i]] = Float64[]
       end
    for file in files
        gene = get_gene(file)
        r = readfunction(joinpath(folder, file))
        writedlm(f,[gene r'],',')
    end
    close(f)
end


function assemble_rates(folder::String,files::Vector,label::String,cond::String,model::String,fish)
    outfile = joinpath(folder,"rates_" * label * "_" * cond * "_" * model * ".csv")
    assemble_files(folder,get_files(files,"rates",label,cond,model),outfile,ratelabels(model,length(split(cond,"-")),fish),readrates)
end

function assemble_measures(folder::String,files,label::String,cond::String,model::String)
    outfile = joinpath(folder,"measures_" * label * "_" * cond * "_" * model *  ".csv")
    header = ["Gene" "Deviance" "LogMaxLikelihood" "WAIC" "SD" "AIC" "Acceptance" "Temperature"]
    assemble_files(folder,get_files(files,"measures",label,cond,model),outfile,header,readmeasures)
end

function assemble_stats(stattype,folder::String,label::String,cond::String,model::String,append::String,header::Bool)
    files = get_files(folder,"param-stats",label,cond,model)
    header = ["Gene"]
    assemble_files(folder,get_files(files,"stats",label,cond,model),outfile,header,readstats)
end


# function assemble_rates(folder::String,label::String,cond::String,model::String,append::String,header::Bool,type,fish)
#     files = get_files(folder,"rates",label,cond,model)
#     println(length(files))
#     outfile = joinpath(folder,"rates_" * label * "_" * cond * "_" * model * append)
#     f = open(outfile,"w")
#     if header
#         r = readrates(joinpath(folder, files[1]),type)
#         writedlm(f,ratelabels(model,length(r),false,fish),',')
#     end
#     for file in files
#         gene = get_gene(file)
#         target = joinpath(folder, file)
#         r = readrates(target,type)
#         writedlm(f,[gene r'],',')
#     end
#     close(f)
# end
#
# function assemble_measures(folder::String,label::String,cond::String,model::String,append::String,header::Bool)
#     files = get_files(folder,"measures",label,cond,model)
#     # label = split(files[1],cond)[1]
#     outfile = joinpath(folder,"measures_" * label * "_" * cond * "_" * model * append)
#     f = open(outfile,"w")
#     if header
#         writedlm(f,["Gene" "Deviance" "LogMaxLikelihood" "WAIC" "SD" "AIC" "Acceptance" "Temperature"],',')
#     end
#     for file in files
#         gene = get_gene(file)
#         target = joinpath(folder, file)
#         r = readmeasures(target)
#         writedlm(f,[gene r],',')
#     end
#     close(f)
# end
# function assemble_stats(stattype,folder::String,label::String,cond::String,model::String,append::String,header::Bool)
#     files = get_files(folder,"param-stats",label,cond,model)
#     # label = split(files[1],cond)[1]
#     outfile = joinpath(folder,"stats_" * label * "_" * cond * "_" * model * append)
#     f = open(outfile,"w")
#     if header
#         # r = readrates(joinpath(folder, files[1]),stattype)
#         writedlm(f,["gene"],',')
#         # writedlm(f,["Gene" "rate01" "sd" "rate10" "sd" "rate12" "sd" "rate21" "sd" "eject" "sd" "yield"],',')
#     end
#     for file in files
#         gene = get_gene(file)
#         target = joinpath(folder, file)
#         r = readstats(target,stattype)
#         println(r)
#         writedlm(f,[gene r'],',')
#     end
#     close(f)
# end

function ratelabels(model,nsets,fish::Bool,sd::Bool=false)
    G = parse(Int,model)
    n = G-1
    Grates = Array{String,2}(undef,1,2*n)
    for i = 0:n-1
        Grates[1,2*i+1] = "Rate$i$(i+1)"
        Grates[1,2*i+2] = "Rate$(i+1)$i"
    end
    lenr = 2*G*nsets + Int(~fish)
    # if fish
    #     lenr =
    #     nsets = div(lenr,2*G)
    # else
    #     nsets = div(lenr-1,2*G)
    # end
    rates = [Grates "Eject" "Decay"]
    for i = 2:nsets
        rates = [rates Grates "Eject" "Decay"]
    end
    if ~fish
        rates = [rates "Yield"]
    end
    if sd
        rates = reshape([rates; fill("SD",1,lenr)],1,2*(lenr))
    end

    return ["Gene" rates]
end

function get_all_rates(file::String,header::Bool)
    r = readdlm(file,',',header=header)
    if header
        r = r[1]
    end
    return r
end



# """
# path_model(type::String,label::String,gene::String,model::String,nalleles,folder,root)
#
# """
# function path_model(type::String,label::String,gene::String,G::Int,R::Int,nalleles::Int,folder,root)
#     file = type * filename(label,gene,G,R,nalleles)
#     joinpath(root,folder,file)
# end
# function path_model(type::String,label::String,gene::String,G::Int,nalleles::Int,folder,root)
#     file = type * filename(label,gene,G,nalleles)
#     joinpath(root,folder,file)
# end

filename(data,model::AbstractGRMmodel) = filename(data.name,data.gene,model.G,model.R,model.nalleles)
filename(data,model::AbstractGMmodel) = filename(data.name,data.gene,model.G,model.nalleles)
filename(label::String,gene::String,G::Int,R::Int,nalleles::Int) = filename(label,gene,"$G"*"$R","$(nalleles)")
filename(label::String,gene,G::Int,nalleles::Int) = filename(label,gene,"$G","$(nalleles)")
filename(label::String,gene::String,model::String,nalleles::String) = "_" * label  * "_" * gene *  "_" * model * "_" * nalleles * txtstr


"""
write_results(file::String,x)
"""
function writeall(path::String,fit,stats,waic,data,temp,model::StochasticGRmodel)
    if ~isdir(path)
        mkpath(path)
    end
    name = filename(data,model)
    write_rates(joinpath(path,"rates" * name ),fit,stats,model)
    write_measures(joinpath(path,"measures" * name),fit,waic,deviance(fit,data,model),temp)
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
function write_measures(file::String,fit::Fit,waic,dev,temp)
    f = open(file,"w")
    writedlm(f,[fit.llml mean(fit.ll) std(fit.ll) quantile(fit.ll,[.025;.5;.975])' waic[1] waic[2] aic(fit)],',')
    writedlm(f,dev,',')
    writedlm(f,[fit.accept fit.total],',')
    writedlm(f,temp,',')
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

function readmeasures(file::String)
    d = readdeviance(file)
    w = readwaic(file)
    a = readaccept(file)
    t = readtemp(file)
    [d[1]; w[1]; w[7]; w[8]; w[9]; a; t[1]]
end

readdeviance(file::String) = readrow(file,2)

readwaic(file::String) = readrow(file,1)

function readaccept(file::String)
    a = readrow(file,3)
    a[1]/a[2]
end

readtemp(file::String) = readrow(file,4)

# function read_stats(file,rows)
#     in = readdlm(file,',')
#     n = sum(in[end,:].!="")
#     in[rows,1:n]
# end

readstats(file::String) = readstats(file,"mean")

function readstats(file::String,stat)
    if stat == "mean"
        m = readmean(file::String)
        return reshape(m,length(m),1)
    else
        return 0
    end
end


function readmean(file::String)
    m = readrow(file,1)
    reshape(m,length(m),1)
end

function readsd(file::String)
    m = readrow(file,2)
    reshape(m,length(m),1)
end

function read_covlogparam(file)
    in = readdlm(file,',')
    n = div((length(in[:,1])-4),4)
    in[end-n+1:end,1:n]
end

function read_covparam(file::String)
    in = readdlm(file,',')
    n = length(in[1,:])
    in[5+2*n:11+2*n,1:n]
end

function read_corparam(file::String)
    in = readdlm(file,',')
    n = length(in[1,:])
    in[5+n:11+n,1:n]
end


"""
readrates(file::String,row::Int)

row
1       maximum likelihood
2       mean
3       median
4       last value of previous run
"""

readrates(file::String) = readrates(file,2)
readrates(file::String,row::Int) = readrow(file,row)

function readrow(file::String,row)
    contents = readdlm(file,',')
    contents[row,:]
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
