# analysis.jl

"""
    make_pca(files::Vector,npcs::Int)

    Compute PCA

"""

function make_combinedpca(files::Vector,npcs::Int)
    m,h = readdlm(files[1],header=true)
    for i in 2:length(files)
        r,h = readdlm(files[i],header=true)
        m = [m r[:,2:end]]
    end
    make_pca(m,npcs)
end

function make_pca(files::Vector,npcs::Int)
    pca = Vector{DataFrame}(undef,length(files))
    for i in 1:length(files)
        pca[i] = make_pca(files[i],npcs)
    end
    return pca
end

function make_pca(file::String,npcs::Int)
    r,h = readdlm(file,header=true)
    make_pca(r,npcs)
end

function make_pca(m::Matrix,npcs::Int)
    M = fit(PCA,float.(m[:,2:end]),maxoutdim = npcs)
    P = projection(M)
    df = DataFrame(Gene = m[:,1])
    i = 1
    for p in eachcol(P)
        insertcols!(df, Symbol("PC$i") => p)
        i += 1
    end
    return df
end

function add_pca(df::DataFrame,files,npcs)
    dpc = make_pca(files,npcs)
    add_pca(df,dpc)
end

add_pca(df::DataFrame,dpc::DataFrame) = leftjoin(df,dpc,on = :Gene)


function isratefile(folder::String)
    files=readdir(folder)
    any(occursin.(".csv",files) .& occursin.("rates",files))
end

function make_dataframe_transient(folder::String,winners::String = "")
    rs = Array{Any,2}(undef,0,8)
    files =readdir(folder)
    n = 0
    for file in files
        if occursin("rate",file)
            t = parse(Float64,split(split(file,"T")[2],"_")[1])
            r,head = readdlm(joinpath(folder,file),',',header=true)
            r = [vcat(r[:,[1,2,3,4,5,10]], r[:,[1,6,7,8,9,10]])  [zeros(size(r)[1]); ones(size(r)[1])]  t*ones(2*size(r)[1])/60.]
            rs = vcat(rs,r)
            n += 1
        end
    end
    if winners != ""
        w = get_winners(winners,2*n)
        return DataFrame(Gene = rs[:,1],on = float.(rs[:,2]),off=float.(rs[:,3]),eject=float.(rs[:,4]),decay=float.(rs[:,5]),yield=float.(rs[:,6]),cond = Int.(rs[:,7]),time = float.(rs[:,8]),winner = w)
    else
        return DataFrame(Gene = rs[:,1],on = float.(rs[:,2]),off=float.(rs[:,3]),eject=float.(rs[:,4]),decay=float.(rs[:,5]),yield=float.(rs[:,6]),cond = Int.(rs[:,7]),time = float.(rs[:,8]))
    end
end


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
"""
    large_deviance(measurefile,threshold)

returns list of genes that have deviance greater than `threshold' in 'measurefile'
"""
large_deviance(measurefile,threshold) = filter_gene(measurefile,"Deviance",threshold)

function deviance(r,cond,n,datafolder,root)
    h,hd = histograms(r,cell,cond,n,datafolder,root)
    deviance(h,hd)
end

function compute_deviance(outfile,ratefile::String,cond,n,datafolder,root)
    f = open(outfile,"w")
    rates,head = readdlm(ratefile,',',header=true)
    for r in eachrow(rates)
        d=deviance(r,cond,n,datafolder,root)
        writedlm(f,[gene d],',')
    end
    close(f)
end


frequency(ron,sd,rdecay) = (ron/rdecay, sd/rdecay)

function burstsize(reject::Float64,roff,covee,covoo,coveo::Float64)
        v = var_ratio(reject,roff,covee,covoo,coveo)
        return reject/roff, sqrt(v)
end

function ratestats(gene,G,folder,cond)
    filestats=joinpath(folder,getfile("param-stats",gene,G,folder,cond)[1])
    filerates = joinpath(folder,getratefile(gene,G,folder,cond)[1])
    rates = readrates(filerates)
    cov = read_covparam(filestats)
    # mu = readmean(filestats)
    # sd = readsd(filestats)
    return rates, cov
end

meanofftime(gene::String,infile,n,method,root) = sum(1 .- offtime(gene,infile,n,method,root))

function meanofftime(r::Vector,n::Int,method::Int)
    if n == 1
        return 1/r[1]
    else
        return sum(1 .- offtime(r,n,method))
    end
end

function offtime(r::Vector,n::Int,method::Int)
    _,_,TI = mat_G_DT(r,n)
    vals,_ = eig_decompose(TI)
    minval = min(minimum(abs.(vals[vals.!=0])),.2)
    offtimeCDF(collect(1.:5/minval),r,n,TI,method)
end

function offtime(gene::String,infile,n,method,root)
    contents,head = readdlm(infile,',',header=true)
    r = float.(contents[gene .== contents[:,1],2:end-1])[1,:]
    offtime(r,n,method)

end

function join_files(file1::String,file2::String,outfile::String,addlabel::Bool=true)
    contents1,head1 = readdlm(file1,',',header=true)   # model G=2
    contents2,head2 = readdlm(file2,',',header=true)   # model G=3
    f = open(outfile,"w")
    if addlabel
        header = vcat(String.(head1[2:end]) .* "_G2",String.(head2[2:end]) .* "_G3")
    else
        header = vcat(String.(head1[2:end]),String.(head2[2:end]))
    end
    header = reshape(permutedims(header),(1,length(head1)+length(head2)-2))
    header = hcat(head1[1],header)
    println(header)
    writedlm(f,header,',')
    for row in 1:size(contents1,1)
        if contents1[row,1] == contents2[row,1]
            contents = hcat(contents1[row:row,2:end],contents2[row:row,2:end])
            contents = reshape(permutedims(contents),(1,size(contents1,2)+size(contents2,2)-2))
            contents = hcat(contents1[row,1],contents)
            writedlm(f,contents,',')
        end
    end
    close(f)
end

function join_files(models::Array,files::Array,outfile::String,addlabel::Bool=true)
    m = length(files)
    contents = Array{Array,1}(undef,m)
    headers = Array{Array,1}(undef,m)
    len = 0
    for i in 1:m
        contents[i],headers[i] = readdlm(files[i],',',header=true)
        len += length(headers[i][2:end])
    end
    f = open(outfile,"w")
    header = Array{String,1}(undef,0)
    for i in 1:m
        if addlabel
            header = vcat(header,String.(headers[i][2:end]) .* "_G$(models[i])")
        else
            header = vcat(header,String.(headers[i][2:end]))
        end
    end
    header = reshape(permutedims(header),(1,len))
    header = hcat(headers[1][1],header)
    println(header)
    writedlm(f,header,',')
    for row in 1:size(contents[1],1)
        content = contents[1][row:row,2:end]
        for i in 1:m-1
            if contents[i][row,1] == contents[i+1][row,1]
                content = hcat(content,contents[i+1][row:row,2:end])
                # content = reshape(permutedims(content),(1,len))
            end
        end
        content = hcat(contents[1][row:row,1],content)
        writedlm(f,[content],',')
    end
    close(f)
end

function sample_non1_genes(infile,n)
    contents,head = readdlm(infile,',',header=true)
    list = Array{String,1}(undef,0)
    for c in eachrow(contents)
        if c[5] != 1
            push!(list,c[1])
        end
    end
    a = StatsBase.sample(list,n,replace=false)
end

function add_best_burst(filein,fileout,filemodel2,filemodel3)
    contents,head = readdlm(filein,',',header=true)
    burst2,head2 = readdlm(filemodel2,',',header=true)
    burst3,head3 = readdlm(filemodel3,',',header=true)
    f = open(fileout,"w")
    head = hcat(head,["mean off period" "bust size"])
    writedlm(f,head,',')
    for row in eachrow(contents)
        if Int(row[end]) == 2
            writedlm(f, hcat(permutedims(row),permutedims(burst2[findfirst(burst2[:,1] .== row[1]),2:3])),',')
        else
            writedlm(f, hcat(permutedims(row),permutedims(burst3[findfirst(burst3[:,1] .== row[1]),2:3])),',')
        end
    end
    close(f)
end

function add_best_occupancy(filein,fileout,filemodel2,filemodel3)
    contents,head = readdlm(filein,',',header=true)
    occupancy2,head2 = readdlm(filemodel2,',',header=true)
    occupancy3,head3 = readdlm(filemodel3,',',header=true)
    f = open(fileout,"w")
    head = hcat(head,["Off -2" "Off -1" "On State" ])
    writedlm(f,head,',')
    for row in eachrow(contents)
        if Int(row[end-2]) == 2
            writedlm(f, hcat(permutedims(row),hcat("NA",permutedims(occupancy2[findfirst(occupancy2[:,1] .== row[1]),2:end]))),',')
        else
            writedlm(f, hcat(permutedims(row),permutedims(occupancy3[findfirst(occupancy3[:,1] .== row[1]),2:end])),',')
        end
    end
    close(f)
end


function prune_file(list,file,outfile,header=true)
    contents,head = readdlm(file,',',header=header)
    f = open(outfile,"w")
    for c in eachrow(contents)
        if c[1] in list
            writedlm(f,[c],',')
        end
    end
    close(f)
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
        gene = get_gene(file1)
        file2 = getratefile(files2,gene)
        outfile = joinpath(outfolder,file2)
        r1 = readrates(joinpath(folder1,file1))
        r2 = readrates(joinpath(folder2,file2))
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
        gene = get_gene(file1)
        file2 = getratefile(files2,gene)
        if file2 != 0
            file2=joinpath(folder2,file2)
        else
            file2=joinpath(folder1,file1)
        end
        name = replace(file1, cond1 => cond1 * "-" * cond2)
        outfile = joinpath(outfolder,name)
        assemble_r(joinpath(folder1,file1),file2,outfile)
    end
end

function  assemble_r(ratefile1,ratefile2,outfile)
    r1 = readrates(ratefile1,2)
    r2 = readrates(ratefile2,2)
    r1[end] = clamp(r1[end],eps(Float64),1-eps(Float64))
    r = vcat(r1[1:end-1],r2[1:end-1],r1[end])
    f = open(outfile,"w")
    writedlm(f,[r],',')
    writedlm(f,[r],',')
    writedlm(f,[r],',')
    writedlm(f,[r],',')
    close(f)
end

function assemble_r(gene,G,folder1,folder2,cond1,cond2,outfolder)
    file1 = getratefile(gene,G,folder1,cond1)[1]
    file2 = getratefile(gene,G,folder2,cond2)[1]
    name = replace(file1, cond1 => cond1 * "-" * cond2)
    println(name)
    outfile = joinpath(outfolder,name)
    println(outfile)
    assemble_r(joinpath(folder1,file1),joinpath(folder2,file2),outfile)
end

function getratefile(files,gene)
    files = files[occursin.("_"*gene*"_",files)]
    if length(files) > 0
        return files[1]
    else
        # println(gene)
        return 0
    end
end

function getratefile(folder,G,cond)
    files = readdir(folder)
    files = files[occursin.("rates_",files)]
    files = files[occursin.("_"*cond*"_",files)]
    files[occursin.("_"*G*"_",files)]
end

getratefile(gene,G,folder,cond) = getfile("rate",gene,G,folder,cond)

function getfile(type,gene::String,G::String,folder,cond)
    files = readdir(folder)
    files = files[occursin.(type,files)]
    files = files[occursin.("_"*gene*"_",files)]
    files = files[occursin.("_"*G*"_",files)]
    files[occursin.("_"*cond*"_",files)]
end

function change_name(folder,oldname,newname)
    files = readdir(folder)
    files = files[occursin.(oldname,files)]
    for file in files
        newfile = replace(file, oldname => newname)
        mv(joinpath(folder,file),joinpath(folder,newfile),force=true)
    end
end

function make_halflife(infile,outfile,col=4)
    f = open(outfile,"w")
    writedlm(f,["Gene" "Halflife"],',')
    contents,rows = readdlm(infile,',',header=true)
    for row = eachrow(contents)
        gene = string(row[1])
        gene = strip(gene,'*')
        h1 = float(row[col])
        h2 = float(row[col+1])
        if h1 > 0 || h2 > 0
            h = (h1 + h2)/(float(h1>0) + float(h2>0))
            writedlm(f,[gene h],',')
        end
    end
    nothing
end

function make_datafiles(infolder,outfolder,label)
    histograms = readdir(infolder)
    if ~isdir(outfolder)
        mkpath(outfolder)
    end
    for file in histograms
        newfile = replace(file,label => "")
        cp(joinpath(infolder,file),joinpath(outfolder,newfile))
    end
end

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
    histograms(r,cell,cond,n::Int,datafolder,root)

"""
function histograms(rin,cell,cond,G::Int,datafolder,fish,root)
    gene = string(rin[1])
    r = float.(rin[2:end])
    data = data_rna(gene,cond,datafolder,fish,"label",root)
    nalleles = alleles(gene,cell,root)
    model = model_rna(r,[],G,nalleles,.01,[],(),fish)
    likelihoodarray(r,data,model)
end

function get_histogram_rna(gene,cond,datafolder,fish,root)
    if fish
        datapath = FISHpath(gene,cond,datafolder,root)
        h = read_fish(datapath,cond,.98)
    else
        datapath = scRNApath(gene,cond,datafolder,root)
        h = read_scrna(datapath,.99)
    end
    normalize_histogram(h)
end

function get_histogram_rna(gene,cond,datafolder,fish)
    if fish
        datapath = FISHpath(gene,cond,datafolder)
        h = read_fish(datapath,cond,.98)
    else
        datapath = scRNApath(gene,cond,datafolder)
        h = read_scrna(datapath,.99)
    end
    normalize_histogram(h)
end


"""
plot_histogram()

functions to plot data and model predicted histograms

"""
function plot_histogram(gene::String,cell::String,G::Int,cond::String,fish::Bool,ratefile::String,datafolder::String,root::String,yield = -1.)
    rates = readdlm(ratefile,',',header=true)
    r = rates[1][findfirst(rates[1][:,1] .== gene)[1],2:end]
    if yield > 0 && ~fish
        r[end] = yield
    end
    data = data_rna(gene,cond,datafolder,fish,"label",root)
    nalleles = alleles(gene,cell,root)
    model = model_rna(r,[],G,nalleles,.01,[],(),fish,0)
    println(typeof(model))
    println(typeof(data))
    m = plot_histogram(data,model)
    return m,data,model
end

function plot_histogram(gene::String,cell::String,G::String,cond::String,fish::Bool,label::String,ratefolder::String,datafolder::String,nsets::Int,root::String,fittedparam = [1],verbose=false)
    data = data_rna(gene,cond,datafolder,fish,label,root)
    model = model_rna(gene,cell,G,fish,.01,fittedparam,(),label,ratefolder,nsets,root,data,verbose)
    m = plot_histogram(data,model)
    return m,data,model
end

function plot_histogram(data::RNAData{Vector{Int64}, Vector{Array}},model::GMlossmodel)
    h=likelihoodarray(model.rates,data,model)
    for i in eachindex(h)
        figure()
        plot(h[i])
        plot(normalize_histogram(data.histRNA[i]))
        savefig(string(i))
    end
    return h
end
function plot_histogram(data::AbstractRNAData{Array{Array,1}},model)
    h=likelihoodarray(model.rates,data,model)
    figure(data.gene)
    for i in eachindex(h)
        plot(h[i])
        plot(normalize_histogram(data.histRNA[i]))
        savefig(string(i))
    end
    return h
end
function plot_histogram(data::AbstractRNAData{Array{Float64,1}},model)
    h=likelihoodfn(get_param(model),data,model)
    plt = plot(h)
    plot!(plt,normalize_histogram(data.histRNA))
    display(plt)
    return h
end

function plot_histogram(data::RNALiveCellData,model::AbstractGRMmodel)
    h=likelihoodtuple(model.rates,data,model)
    plt1 = plot(h[1])
    plot!(plt1,normalize_histogram(data.OFF))
    plt2 = plot(h[2])
    plot!(plt2,normalize_histogram(data.ON))
    plt3 = plot(h[3])
    plot!(plt3,normalize_histogram(data.histRNA))
    plt = plot(plt1,plt2,plt3,layout = (3,1))
    display(plt)
    return h
end

function plot_histogram(data::TransientRNAData,model::AbstractGMmodel)
    h=likelihoodarray(model.rates,data,model)
    for i in eachindex(h)
        figure(data.gene *":T" * "$(data.time[i])")
        plot(h[i])
        plot(normalize_histogram(data.histRNA[i]))
    end
    return h
end

function plot_histogram(data::RNAData{T1,T2},model::AbstractGMmodel,save = false) where {T1 <: Array, T2 <: Array}
    m=likelihoodarray(model.rates,data,model)
    println("*")
    for i in eachindex(m)
        plt=plot(m[i])
        plot!(normalize_histogram(data.histRNA[i]),show=true)
        if save
            savefig()
        else
            display(plt)
        end
    end
    println(deviance(data,model))
    println(loglikelihood(get_param(model),data,model)[1])
    return m
end

function plot_histogram(data::RNAData,model::AbstractGMmodel)
    h=likelihoodfn(get_param(model),data,model)
    plt=plot(h)
    plot!(normalize_histogram(data.histRNA))
    display(plt)
    return h
end

function plot_model(r,n,nhist,nalleles,yield)
    h= steady_state(r[1:2*n+2],yield,n,nhist,nalleles)
    plt = plot(h)
    display(plt)
    return h
end

function plot_model(r,n,nhist,nalleles)
    h= steady_state(r[1:2*n+2],n,nhist,nalleles)
    plt = plot(h)
    display(plt)
    return h
end
