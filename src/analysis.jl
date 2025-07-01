# This file is part of StochasticGene.jl   

# analysis.jl

"""
    make_dataframes(resultfolder::String,datapath::String, assemble=true, multicond=false, datatype="rna")

"""
function make_dataframes(resultfolder::String, datapath::String, assemble=true, multicond=false, datatype="rna")
    if assemble
        assemble_all(resultfolder, multicond)
    end
    files = get_ratesummaryfiles(resultfolder)
    parts = fields.(files)
    models = get_models(parts)
    labels = get_labels(parts)
    df = Vector{Vector}(undef, 0)
    for label in labels
        lfiles = files[label.==get_label.(files)]
        dfl = Vector{Tuple}(undef, 0)
        for model in models
            if length(parse_model(model)) > 1
                println("Not G model")
            else
                mfiles = lfiles[model.==get_model.(lfiles)]
                dfm = Vector{DataFrame}(undef, 0)
                for i in eachindex(mfiles)
                    push!(dfm, make_dataframe(joinpath(resultfolder, mfiles[i]), datapath, multicond, datatype))
                end
                push!(dfl, ("Summary_$(label)_$(model).csv", stack_dataframe(dfm)))
            end
        end
        push!(df, dfl)
    end
    return df
end

"""
    statfile_from_ratefile(ratefile)

TBW
"""
statfile_from_ratefile(ratefile) = replace(ratefile, "rates_" => "stats_")

"""
    make_dataframe(ratefile::String, datapath::String)

TBW
"""
function make_dataframe(ratefile::String, datapath::String, multicond::Bool, datatype::String)
    df = read_dataframe(ratefile)
    df2 = read_dataframe(statfile_from_ratefile(ratefile))
    df = leftjoin(df, df2, on=:Gene, makeunique=true)
    filename = splitpath(ratefile)[end]
    parts = fields(filename)
    model = parse_model(parts.model)
    if length(model) > 1
        println("Not G model")
    else
        G = parse(Int, parts.model)
        insertcols!(df, :Model => fill(G, size(df, 1)))
        df = stack_dataframe(df, G, parts.cond, multicond)
        add_moments!(df, datapath, datatype)
    end
end

"""
    parse_model(name::String)

TBW
"""
function parse_model(name::String)
    d = parse(Int, name)
    d > 9 && (d = digits(d))
    d
end

"""
    augment_dataframe(df, resultfolder)

TBW
"""
function augment_dataframe(df, resultfolder)
    dfc = copy(df)
    G = dfc[1, :Model]
    if G > 1
        dfc = add_burstsize(dfc, resultfolder)
    end
    if G == 2
        add_modelmoments!(dfc)
    end
    dfc = add_measures(dfc, resultfolder, string(G))
    add_residenceprob!(dfc)
    dfc
end

"""
    make_measure_df(resultfolder::String, G::String)

TBW
"""
function make_measure_df(resultfolder::String, G::String)
    files = get_measuresummaryfiles(resultfolder)
    df = Vector{DataFrame}(undef, 0)
    for file in files
        parts = fields(file)
        if parts.model == G
            df0 = read_dataframe(joinpath(resultfolder, file))
            if ~isempty(parts.cond)
                insertcols!(df0, :Condition => parts.cond)
            end
            push!(df, df0)
        end
    end
    stack_dataframe(df)
end

"""
    join_cols(d)

TBW
"""
function join_cols(d)
    if ismissing(d.Condition[1])
        return [:Gene]
    else
        return [:Gene, :Condition]
    end
end

"""
    add_measures(df, resultfolder::String, G)

TBW
"""
function add_measures(df, resultfolder::String, G)
    dm = make_measure_df(resultfolder, G)
    leftjoin(df, dm, on=join_cols(df))
end

"""
    add_mean!(df::DataFrame, datapath)

TBW
"""
function add_mean!(df::DataFrame, datapath)
    # root = string(split(abspath(datapath),"data")[1])
    m = Vector{Float64}(undef, length(df.Gene))
    i = 1
    for gene in df.Gene
        m[i] = mean_histogram(get_histogram_rna(string(gene), df[i, :Condition], datapath))
        i += 1
    end
    insertcols!(df, :Expression => m)
end

"""
    add_moments!(df::DataFrame, datapath)

TBW
"""
function add_moments!(df::DataFrame, datapath, datatype::String)
    m = Vector{Float64}(undef, length(df.Gene))
    v = similar(m)
    t = similar(m)
    i = 1
    if datatype == "rna"
        for gene in df.Gene
            _, h = read_rna(gene, df[i, :Condition], datapath)
            m[i] = mean_histogram(h)
            v[i] = var_histogram(h)
            t[i] = moment_histogram(h, 3)
            i += 1
        end
        insertcols!(df, :Expression => m, :Variance => v, :ThirdMoment => t)
    elseif datatype == "rnacount"
        for gene in df.Gene
            c, _, _ = read_rnacount(gene, df[i, :Condition], datapath)
            m[i] = mean(c)
            v[i] = var(c)
            t[i] = moment(c, 3)
            i += 1
        end
        insertcols!(df, :Expression => m, :Variance => v, :ThirdMoment => t)
    end
end

"""
    add_modelmoments!(df::DataFrame)

TBW
"""
function add_modelmoments!(df::DataFrame)
    m = Vector{Float64}(undef, length(df.Gene))
    v = similar(m)
    i = 1
    for gene in df.Gene
        m[i] = model2_mean(df.Rate12[i], df.Rate21[i], df.Eject[i], df.Decay[i], df.Nalleles[i])
        v[i] = model2_variance(df.Rate12[i], df.Rate21[i], df.Eject[i], df.Decay[i], df.Nalleles[i])
        i += 1
    end
    insertcols!(df, :Model_Expression => m, :Model_Variance => v)
end

"""
    add_residenceprob!(df::DataFrame)

add residence probability of G states to dataframe
Convention changed from original code. G states are now labeled starting from 1 rather than zero
"""
function add_residenceprob!(df::DataFrame)
    n = df.Model[1] - 1
    N = length(df.Gene)
    g = Matrix{Float64}(undef, N, n + 1)
    r = Vector{Float64}(undef, 2 * n)
    for i in 1:N
        for j in 1:n
            Symbol("Rate$j$(j-1)")
            Symbol("Rate$(j-1)$j")
            r[2*j-1] = df[i, Symbol("Rate$j$(j-1)")]
            r[2*j] = df[i, Symbol("Rate$(j-1)$j")]
            g[i, :] = residenceprob_G(r, n + 1)
        end
    end
    for j in 1:n+1
        insertcols!(df, Symbol("ProbG$(j-1)") => g[:, j])
    end
end

"""
    add_burstsize(df, resultfolder::String, cols::Vector{Symbol}=join_cols(df))

TBW
"""
add_burstsize(df, resultfolder::String, cols::Vector{Symbol}=join_cols(df)) = add_burstsize(df, make_burst_df(resultfolder), cols)

"""
    add_burstsize(df, db, cols)

TBW
"""
function add_burstsize(df, db, cols)
    if ismissing(df[1, :Condition])
        leftjoin(df, db[:, [:BurstMean, :BurstSD, :BurstMedian, :BurstMAD, :Gene]], on=cols)
    else
        leftjoin(df, db[:, [:BurstMean, :BurstSD, :BurstMedian, :BurstMAD, :Gene, :Condition]], on=cols)
    end
end

"""
    make_burst_df(resultfolder::String)

TBW
"""
function make_burst_df(resultfolder::String)
    files = get_burstsummaryfiles(resultfolder)
    df = Vector{DataFrame}(undef, 0)
    for file in files
        parts = fields(file)
        df0 = read_dataframe(joinpath(resultfolder, file))
        if ~isempty(parts.cond)
            insertcols!(df0, :Condition => parts.cond)
        end
        push!(df, df0)
    end
    stack_dataframe(df)
end


# add_time(csvfile::String,timestamp) = CSV.write(csvfile,add_time!(read_dataframe(csvfile),timestamp))

"""
    add_time!(df::DataFrame, timestamp)

TBW
"""
add_time!(df::DataFrame, timestamp) = insertcols!(df, :Time => timestamp)

"""
    stack_dataframe(df,G,cond)
    stack_dataframe(df2::Vector{DataFrame})

"""
stack_dataframe(df, G, cond, multicond) = stack_dataframe(separate_dataframe(df, G, cond, multicond))

function stack_dataframe(df2::Vector{DataFrame})
    df = df2[1]
    for i in 2:length(df2)
        df = vcat(df, df2[i])
    end
    return df
end

"""
    separate_dataframe(df, G, cond)

TBW
"""
function separate_dataframe(df, G, cond, multicond)
    println("cond: ", cond)
    conds = split_conditions(cond, multicond)
    nsets = length(conds)
    df2 = Vector{DataFrame}(undef, nsets)
    for i in 1:nsets
        df2[i] = df[:, [1; 2*G*(i-1)+2:2*G*i+1; 2*G*nsets+2:end]]
        if multicond
            rename!(x -> split(x, "_")[1], df2[i])
        end
        insertcols!(df2[i], :Condition => fill(string(conds[i]), size(df, 1)))
    end
    return df2
end

function separate_dataframe(df, G)

    df2 = Vector{DataFrame}(undef, nsets)
    for i in 1:nsets
        df2[i] = df[:, [1; 2*G*(i-1)+2:2*G*i+1; 2*G*nsets+2:end]]
    end
    return df2
end

# function residuals(df)
#     R01 = lm(@formula(log(Rate01) ~  Time*log(Decay)+PC1+PC2+PC3+PC4+PC5),df[(df.Winner .> 1) .& (df.Expression .> 0.),:])
#
# end

"""
    make_Zscore_dataframe(df, conditions::Vector)

TBW
"""
function make_Zscore_dataframe(df, conditions::Vector)
    dfc = Array{DataFrame}(undef, length(conditions))
    for i in eachindex(dfc)
        dfc[i] = df[df.Condition.==conditions[i], :]
    end
    df = leftjoin(dfc[1], dfc[2], on=[:Gene, :Time], makeunique=true)
    add_Zscore!(df)
    add_MomentChange!(df)
    df[:, [:Gene, :Time, :ZRate01, :ZRate10, :ZEject, :dExpression, :dVariance, :dThirdMoment, :Winner]]
end

"""
    add_MomentChange!(df)

TBW
"""
function add_MomentChange!(df)
    insertcols!(df, :dExpression => df.Expression_1 .- df.Expression)
    insertcols!(df, :dVariance => df.Variance_1 .- df.Variance)
    insertcols!(df, :dThirdMoment => df.ThirdMoment_1 .- df.ThirdMoment)
end

"""
    add_Zscore!(df)

TBW
"""
function add_Zscore!(df)
    insertcols!(df, :ZRate01 => Difference_Zscore.(df.Rate01_1, df.Rate01, df.Rate01SD_1, df.Rate01SD))
    insertcols!(df, :ZRate10 => Difference_Zscore.(df.Rate10_1, df.Rate10, df.Rate10SD_1, df.Rate10SD))
    insertcols!(df, :ZEject => Difference_Zscore.(df.Eject_1, df.Eject, df.EjectSD, df.EjectSD_1))
end

function add_Zscore_class!(df, threshold=2.0)
    insertcols!(df, :On => classify_Zscore.(df.ZRate01, threshold))
    insertcols!(df, :Off => classify_Zscore.(df.ZRate10, threshold))
    insertcols!(df, :Eject => classify_Zscore.(df.ZEject, threshold))
    insertcols!(df, :OnOffEject => df.On .* df.Off .* df.Eject)
end

function classify_Zscore(Z, threshold)
    if Z .> threshold
        return "Up"
    elseif Z .< -threshold
        return "Down"
    else
        return "Null"
    end
end

"""
    best_AIC(folder::String)

TBW
"""
best_AIC(folder::String) = best_measure(folder, :AIC)

"""
    best_WAIC(folder::String)

TBW
"""
best_WAIC(folder::String) = best_measure(folder, :WAIC)

"""
    best_measure(folder::String, measure::Symbol)

TBW
"""
function best_measure(folder::String, measure::Symbol)
    files = get_measuresummaryfiles(folder)
    parts = fields.(files)
    labels = get_labels(parts)
    conds = get_conds(parts)
    df = Vector{Tuple{String,DataFrame}}(undef, 0)
    for label in labels
        for cond in conds
            lfiles = files[(label.==get_label.(files)).&(cond.==get_cond.(files))]
            dm = Vector{DataFrame}(undef, length(lfiles))
            for i in eachindex(lfiles)
                if isfile(joinpath(folder, lfiles[i]))
                    dm[i] = read_dataframe(joinpath(folder, lfiles[i]))
                    insertcols!(dm[i], :Model => fill(get_model(lfiles[i]), size(dm[i], 1)))
                end
            end
            # println(best_measure(dm,measure))
            bm = best_measure(dm, measure)
            ~isempty(bm) && push!(df, ("Winners_$(label)_$(cond)_$(string(measure)).csv", bm))
        end
    end
    return df
end

function best_measure(dfs::Vector, measure::Symbol)
    df = DataFrame(:Gene => [], :Winner => [], measure => [])
    ngenes = Int[]
    for d in dfs
        ngenes = push!(ngenes, length(d[:, :Gene]))
    end
    if ~isempty(ngenes)
        others = setdiff(eachindex(dfs), argmax(ngenes))
        for d in eachrow(dfs[argmax(ngenes)])
            l = d[measure]
            model = d[:Model]
            for k in others
                dc = dfs[k][dfs[k].Gene.==d.Gene, :]
                if ~isempty(dc)
                    if dc[1, measure] < l
                        l = dc[1, measure]
                        model = dc[1, :Model]
                    end
                end
            end
            push!(df, Dict(:Gene => d.Gene, :Winner => model, measure => l))
        end
    end
    return df
end


"""
    add_pca(df::DataFrame,files::Vector,npcs::Int,conds::Vector)

    make PCAs from files and add to dataframe df

"""

add_pca(df::DataFrame, files::Vector, npcs::Int, conds::Vector) = add_pca(df, make_pca(files, npcs, conds), makeunique=true)

add_pca(df::DataFrame, dpc::DataFrame) = leftjoin(df, dpc, on=[:Condition, :Gene, :Time], makeunique=true)


"""
    make_pca(files::Vector,npcs::Int)

    Compute PCA

"""

function make_pca(files::Vector, npcs::Int, conds::Vector, times::Vector)
    pca = make_pca(files[1], npcs, conds[1], times[1])
    for i in 2:length(files)
        pca = [pca; make_pca(files[i], npcs, conds[i], times[i])]
    end
    return pca
end

function make_pca(files::Vector, npcs::Int, conds::Vector)
    pca = make_pca(files[1], npcs, conds[1])
    for i in 2:length(files)
        pca = [pca; make_pca(files[i], npcs, conds[i])]
    end
    return pca
end

function make_pca(file::String, npcs::Int, cond)
    r, h = readdlm(file, header=true)
    df = make_pca(r, npcs)
    insertcols!(df, :Condition => cond)
end

function make_pca(file::String, npcs::Int, cond, time)
    r, h = readdlm(file, header=true)
    df = make_pca(r, npcs)
    insertcols!(df, :Condition => cond, :Time => time)
end

function make_pca(m::Matrix, npcs::Int)
    M = MultivariateStats.fit(PCA, float.(m[:, 2:end]), maxoutdim=npcs)
    make_pca(M, m[:, 1])
end

function make_pca(M::PCA, genes)
    P = projection(M)
    df = DataFrame(Gene=genes)
    i = 1
    for p in eachcol(P)
        insertcols!(df, Symbol("PC$i") => p)
        i += 1
    end
    return df
end

function make_combinedpca(files::Vector, npcs::Int)
    m, h = readdlm(files[1], header=true)
    for i in 2:length(files)
        r, h = readdlm(files[i], header=true)
        m = [m r[:, 2:end]]
    end
    make_pca(m, npcs)
end


function make_dataframe_transient(folder::String, winners::String="")
    rs = Array{Any,2}(undef, 0, 8)
    files = readdir(folder)
    n = 0
    for file in files
        if occursin("rate", file)
            t = parse(Float64, split(split(file, "T")[2], "_")[1])
            r, head = readdlm(joinpath(folder, file), ',', header=true)
            r = [vcat(r[:, [1, 2, 3, 4, 5, 10]], r[:, [1, 6, 7, 8, 9, 10]]) [zeros(size(r)[1]); ones(size(r)[1])] t * ones(2 * size(r)[1]) / 60.0]
            rs = vcat(rs, r)
            n += 1
        end
    end
    if winners != ""
        w = get_winners(winners, 2 * n)
        return DataFrame(Gene=rs[:, 1], on=float.(rs[:, 2]), off=float.(rs[:, 3]), eject=float.(rs[:, 4]), decay=float.(rs[:, 5]), yield=float.(rs[:, 6]), cond=Int.(rs[:, 7]), time=float.(rs[:, 8]), winner=w)
    else
        return DataFrame(Gene=rs[:, 1], on=float.(rs[:, 2]), off=float.(rs[:, 3]), eject=float.(rs[:, 4]), decay=float.(rs[:, 5]), yield=float.(rs[:, 6]), cond=Int.(rs[:, 7]), time=float.(rs[:, 8]))
    end
end

function filter_gene(measurefile, measure, threshold)
    genes = Vector{String}(undef, 0)
    measures, header = readdlm(measurefile, ',', header=true)
    println("length of measures: ", length(measures[:, 1]))
    col = findfirst(header[1, :] .== measure)
    for d in eachrow(measures)
        if d[col] > threshold || isnan(d[col])
            push!(genes, d[1])
        end
    end
    println("length of genes: ", length(genes))
    return genes
end

function filter_gene_nan(measurefile, measure)
    genes = Vector{String}(undef, 0)
    measures, header = readdlm(measurefile, ',', header=true)
    println("length of measures: ", length(measures[:, 1]))
    col = findfirst(header[1, :] .== measure)
    for d in eachrow(measures)
        if isnan(d[col])
            push!(genes, d[1])
        end
    end
    println("length of genes: ", length(genes))
    return genes
end
"""
large_rhat(measureile,threshold)

"""
large_rhat(measurefile, threshold) = filter_gene(measurefile, "Rhat", threshold)
"""
    large_deviance(measurefile,threshold)

returns list of genes that have deviance greater than `threshold' in 'measurefile'
"""
large_deviance(measurefile, threshold) = filter_gene(measurefile, "Deviance", threshold)

function deviance(r, cond, n, datapath, root)
    h, hd = histograms(r, cell, cond, n, datapath, root)
    deviance(h, hd)
end

function compute_deviance(outfile, ratefile::String, cond, n, datapath, root)
    f = open(outfile, "w")
    rates, head = readdlm(ratefile, ',', header=true)
    for r in eachrow(rates)
        d = deviance(r, cond, n, datapath, root)
        writedlm(f, [gene d], ',')
    end
    close(f)
end


"""
    deviance(fits, data, model)

return deviance
"""
deviance(fits, data, model) = -1.0

function deviance(fits, data::AbstractHistogramData, model)
    predictions = predictedfn(fits.parml, data, model)
    deviance(log.(max.(predictions, eps())), datapdf(data))
end

deviance(fits, data::RNACountData, model) = -1


"""
    deviance(fits, data::AbstractTraceData, model)

return max ll normalized by number of trace frames
"""
function deviance(fits, data::AbstractTraceData, model)
    if isempty(data.trace[1])
        return -1.0
    else
        return fits.llml / sum(sum.(data.trace[1]))
    end
end

function deviance(fits, data::AbstractTraceData, model::AbstractGRSMmodel{TraitType}) where {TraitType<:NamedTuple}
    if hastrait(model, :coupling)
        # Coupling-specific logic
        s = 0
        for t in data.trace[1]
            s += sum(sum.(t))
        end
        return fits.llml / s
    else
        # Default logic
        s = 0
        for t in data.trace[1]
            s += sum(sum.(t))
        end
        return fits.llml / s
    end
end



"""
    deviance(data::AbstractHistogramData, model::AbstractGeneTransitionModel)


"""
function deviance(data::AbstractHistogramData, model::AbstractGeneTransitionModel)
    h = predictedfn(model.rates[model.fittedparam], data, model)
    println("h: ", h)
    deviance(h, datapdf(data))
end

"""
    deviance(logpredictions::Array,data::AbstractHistogramData)

return deviance

use log of data histogram as loglikelihood of overfit model
"""
deviance(logpredictions::Array, data::AbstractHistogramData) = deviance(logpredictions, datapdf(data))

deviance(logpredictions::Array, hist::Array) = 2 * hist' * (log.(max.(hist, eps())) - logpredictions)



"""
    frequency(ron, sd, rdecay)

return on rate / decay rate, sd
"""
frequency(ron, sd, rdecay) = (ron / rdecay, sd / rdecay)


"""
    burstsize(reject::Float64, roff, covee, covoo, coveo::Float64)

return eject rate / off rate, sd
"""
function burstsize(reject::Float64, roff, covee, covoo, coveo::Float64)
    v = var_ratio(reject, roff, covee, covoo, coveo)
    return reject / roff, sqrt(v)
end

function ratestats(gene, G, folder, cond)
    filestats = joinpath(folder, getfile("param-stats", gene, G, folder, cond)[1])
    filerates = joinpath(folder, getratefile(gene, G, folder, cond)[1])
    rates = readrates(filerates)
    cov = read_covparam(filestats)
    # mu = readmean(filestats)
    # sd = readsd(filestats)
    return rates, cov
end

meanofftime(gene::String, infile, n, method, root) = sum(1 .- offtime(gene, infile, n, method, root))

function meanofftime(r::Vector, n::Int, method)
    if n == 1
        return 1 / r[1]
    else
        return sum(1 .- offtime(r, n, method))
    end
end

function offtime(r::Vector, n::Int, method)
    _, _, TI = mat_G_DT(r, n)
    vals, _ = eig_decompose(TI)
    minval = min(minimum(abs.(vals[vals.!=0])), 0.2)
    offtimeCDF(collect(1.0:5/minval), r, n, TI, method)
end

function offtime(gene::String, infile, n, method, root)
    contents, head = readdlm(infile, ',', header=true)
    r = float.(contents[gene.==contents[:, 1], 2:end-1])[1, :]
    offtime(r, n, method)

end

function join_files(file1::String, file2::String, outfile::String, addlabel::Bool=true)
    contents1, head1 = readdlm(file1, ',', header=true)   # model G=2
    contents2, head2 = readdlm(file2, ',', header=true)   # model G=3
    f = open(outfile, "w")
    if addlabel
        header = vcat(String.(head1[2:end]) .* "_G2", String.(head2[2:end]) .* "_G3")
    else
        header = vcat(String.(head1[2:end]), String.(head2[2:end]))
    end
    header = reshape(permutedims(header), (1, length(head1) + length(head2) - 2))
    header = hcat(head1[1], header)
    println("header: ", header)
    writedlm(f, header, ',')
    for row in 1:size(contents1, 1)
        if contents1[row, 1] == contents2[row, 1]
            contents = hcat(contents1[row:row, 2:end], contents2[row:row, 2:end])
            contents = reshape(permutedims(contents), (1, size(contents1, 2) + size(contents2, 2) - 2))
            contents = hcat(contents1[row, 1], contents)
            writedlm(f, contents, ',')
        end
    end
    close(f)
end

function join_files(models::Array, files::Array, outfile::String, addlabel::Bool=true)
    m = length(files)
    contents = Array{Array,1}(undef, m)
    headers = Array{Array,1}(undef, m)
    len = 0
    for i in 1:m
        contents[i], headers[i] = readdlm(files[i], ',', header=true)
        len += length(headers[i][2:end])
    end
    f = open(outfile, "w")
    header = Array{String,1}(undef, 0)
    for i in 1:m
        if addlabel
            header = vcat(header, String.(headers[i][2:end]) .* "_G$(models[i])")
        else
            header = vcat(header, String.(headers[i][2:end]))
        end
    end
    header = reshape(permutedims(header), (1, len))
    header = hcat(headers[1][1], header)
    println("header: ", header)
    writedlm(f, header, ',')
    for row in 1:size(contents[1], 1)
        content = contents[1][row:row, 2:end]
        for i in 1:m-1
            if contents[i][row, 1] == contents[i+1][row, 1]
                content = hcat(content, contents[i+1][row:row, 2:end])
                # content = reshape(permutedims(content),(1,len))
            end
        end
        content = hcat(contents[1][row:row, 1], content)
        writedlm(f, [content], ',')
    end
    close(f)
end

function sample_non1_genes(infile, n)
    contents, head = readdlm(infile, ',', header=true)
    list = Array{String,1}(undef, 0)
    for c in eachrow(contents)
        if c[5] != 1
            push!(list, c[1])
        end
    end
    a = StatsBase.sample(list, n, replace=false)
end

function add_best_burst(filein, fileout, filemodel2, filemodel3)
    contents, head = readdlm(filein, ',', header=true)
    burst2, head2 = readdlm(filemodel2, ',', header=true)
    burst3, head3 = readdlm(filemodel3, ',', header=true)
    f = open(fileout, "w")
    head = hcat(head, ["mean off period" "bust size"])
    writedlm(f, head, ',')
    for row in eachrow(contents)
        if Int(row[end]) == 2
            writedlm(f, hcat(permutedims(row), permutedims(burst2[findfirst(burst2[:, 1] .== row[1]), 2:3])), ',')
        else
            writedlm(f, hcat(permutedims(row), permutedims(burst3[findfirst(burst3[:, 1] .== row[1]), 2:3])), ',')
        end
    end
    close(f)
end

function add_best_occupancy(filein, fileout, filemodel2, filemodel3)
    contents, head = readdlm(filein, ',', header=true)
    occupancy2, head2 = readdlm(filemodel2, ',', header=true)
    occupancy3, head3 = readdlm(filemodel3, ',', header=true)
    f = open(fileout, "w")
    head = hcat(head, ["Off -2" "Off -1" "On State"])
    writedlm(f, head, ',')
    for row in eachrow(contents)
        if Int(row[end-2]) == 2
            writedlm(f, hcat(permutedims(row), hcat("NA", permutedims(occupancy2[findfirst(occupancy2[:, 1] .== row[1]), 2:end]))), ',')
        else
            writedlm(f, hcat(permutedims(row), permutedims(occupancy3[findfirst(occupancy3[:, 1] .== row[1]), 2:end])), ',')
        end
    end
    close(f)
end

function prune_file(list, file, outfile, header=true)
    contents, head = readdlm(file, ',', header=header)
    f = open(outfile, "w")
    for c in eachrow(contents)
        if c[1] in list
            writedlm(f, [c], ',')
        end
    end
    close(f)
end

function replace_yield(G, folder1, folder2, cond1, cond2, outfolder)
    if typeof(G) <: Number
        G = string(G)
    end
    if ~isdir(outfolder)
        mkpath(outfolder)
    end
    files1 = getratefile(folder1, G, cond1)
    files2 = getratefile(folder2, G, cond2)
    for file1 in files1
        gene = get_gene(file1)
        file2 = getratefile(files2, gene)
        outfile = joinpath(outfolder, file2)
        r1 = readrates(joinpath(folder1, file1))
        r2 = readrates(joinpath(folder2, file2))
        r2[end] = r1[end]
        f = open(outfile, "w")
        writedlm(f, [r2], ',')
        close(f)
    end

end


"""
assemble_r(G,folder1,folder2,cond1,cond2,outfolder)

Combine rates from two separate fits into a single rate vector

"""

function assemble_r(G, folder1, folder2, cond1, cond2, outfolder)
    if typeof(G) <: Number
        G = string(G)
    end
    if ~isdir(outfolder)
        mkpath(outfolder)
    end
    files1 = getratefile(folder1, G, cond1)
    files2 = getratefile(folder2, G, cond2)
    for file1 in files1
        gene = get_gene(file1)
        file2 = getratefile(files2, gene)
        if file2 != 0
            file2 = joinpath(folder2, file2)
        else
            file2 = joinpath(folder1, file1)
        end
        name = replace(file1, cond1 => cond1 * "-" * cond2)
        outfile = joinpath(outfolder, name)
        assemble_r(joinpath(folder1, file1), file2, outfile)
    end
end

function assemble_r(ratefile1, ratefile2, outfile)
    r1 = readrates(ratefile1, 2)
    r2 = readrates(ratefile2, 2)
    r1[end] = clamp(r1[end], eps(Float64), 1 - eps(Float64))
    r = vcat(r1[1:end-1], r2[1:end-1], r1[end])
    f = open(outfile, "w")
    writedlm(f, [r], ',')
    writedlm(f, [r], ',')
    writedlm(f, [r], ',')
    writedlm(f, [r], ',')
    close(f)
end

function assemble_r(gene, G, folder1, folder2, cond1, cond2, outfolder)
    file1 = getratefile(gene, G, folder1, cond1)[1]
    file2 = getratefile(gene, G, folder2, cond2)[1]
    name = replace(file1, cond1 => cond1 * "-" * cond2)
    println("name: ", name)
    outfile = joinpath(outfolder, name)
    println("outfile: ", outfile)
    assemble_r(joinpath(folder1, file1), joinpath(folder2, file2), outfile)
end

function getratefile(files, gene)
    files = files[occursin.("_" * gene * "_", files)]
    if length(files) > 0
        return files[1]
    else
        # println(gene)
        return 0
    end
end

function getratefile(folder, G, cond)
    files = readdir(folder)
    files = files[occursin.("rates_", files)]
    files = files[occursin.("_" * cond * "_", files)]
    files[occursin.("_" * G * "_", files)]
end

getratefile(gene, G, folder, cond) = getfile("rate", gene, G, folder, cond)

function getfile(type, gene::String, G::String, folder, cond)
    files = readdir(folder)
    files = files[occursin.(type, files)]
    files = files[occursin.("_" * gene * "_", files)]
    files = files[occursin.("_" * G * "_", files)]
    files[occursin.("_" * cond * "_", files)]
end

function change_name(folder, oldname, newname)
    files = readdir(folder)
    files = files[occursin.(oldname, files)]
    for file in files
        newfile = replace(file, oldname => newname)
        mv(joinpath(folder, file), joinpath(folder, newfile), force=true)
    end
end

function make_halflife(infile, outfile, col=4)
    f = open(outfile, "w")
    writedlm(f, ["Gene" "Halflife"], ',')
    contents, rows = readdlm(infile, ',', header=true)
    for row = eachrow(contents)
        gene = string(row[1])
        gene = strip(gene, '*')
        h1 = float(row[col])
        h2 = float(row[col+1])
        if h1 > 0 || h2 > 0
            h = (h1 + h2) / (float(h1 > 0) + float(h2 > 0))
            writedlm(f, [gene h], ',')
        end
    end
    nothing
end

function make_datafiles(infolder, outfolder, label)
    histograms = readdir(infolder)
    if ~isdir(outfolder)
        mkpath(outfolder)
    end
    for file in histograms
        newfile = replace(file, label => "")
        cp(joinpath(infolder, file), joinpath(outfolder, newfile))
    end
end


"""
    histograms(r,cell,cond,n::Int,datapath,root)

"""
function histograms(rin, cell, cond, G::Int, datapath, root)
    fish = false
    gene = string(rin[1])
    r = float.(rin[2:end])
    data = data_rna(gene, cond, datapath, fish, "label", root)
    nalleles = alleles(gene, cell, root)
    model = model_rna(r, [], G, nalleles, 0.01, [], (), fish)
    predictedarray(r, data, model)
end

function get_histogram_rna(gene, datacond, datapath)
    _, h = read_rna(gene, datacond, datapath)
    normalize_histogram(h)
end



"""
    write_ONOFFhistograms(r, transitions, G, R, S, insertstep, bins; outfile::String="")

simulations and master equation solutions of dwell time histograms
"""
function write_ONOFFhistograms(r, transitions, G, R, S, insertstep, bins; outfile::String="", simulate=false)
    onstates = on_states(G, R, S, insertstep)
    components = TAIComponents(transitions, G, R, S, insertstep, onstates, "")
    T = make_mat_T(components, r)
    TA = make_mat_TA(components, r)
    TI = make_mat_TI(components, r)
    OFF, ON = offonPDF(bins, r, T, TA, TI, components.nT, components.elementsT, onstates)
    if simulate
        hs = simulator(r, transitions, G, R, S, insertstep, bins=bins)
        df = DataFrame(time=bins, ON=ON, OFF=OFF, SimON=hs[2] / sum(hs[2]), SimOFF=hs[3] / sum(hs[3]))
    else
        df = DataFrame(time=bins, ON=ON, OFF=OFF)
    end
    if ~isempty(outfile)
        CSV.write(outfile, df)
    end
    return df
end

function write_ONOFFhistograms(folder::String, bins=collect(1.0:200.0))
    files = get_resultfiles(folder)
    for (root, dirs, files) in walkdir(folder)
        for f in files
            if occursin("rates", f)
                println(joinpath(root, f))
                parts = fields(f)
                G, R, S, insertstep = decompose_model(parts.model)
                r = readrates(joinpath(root, f))
                out = joinpath(root, replace(f, "rates" => "ONOFF", ".txt" => ".csv"))
                transitions = get_transitions(G, parts.label)
                write_ONOFFhistograms(r, transitions, G, R, S, insertstep, bins, outfile=out)
            end
        end
    end
end


function write_RNAhistogram(r, transitions, G::Int, R::Int, decayrate::Float64, nalleles::Int, nRNA::Int; outfile::String="", splicetype="", ejectnumber=1)
    mcomponents = MComponents(transitions, G, R, nRNA, decayrate, splicetype, ejectnumber)
    df = DataFrame(Freq=predictedRNA(r, mcomponents, nalleles, nRNA), Bin=collect(0:nRNA-1))
    if ~isempty(outfile)
        CSV.write(outfile, df)
    end
    return df
end

function write_RNAhistogram(folder, nRNA; ejectnumber=1)
    files = get_resultfiles(folder)
    for (root, dirs, files) in walkdir(folder)
        for f in files
            if occursin("rates", f) && !occursin("joint", f)
                println(f)
                parts = fields(f)
                G, R, S, insertstep = decompose_model(parts.model)
                nalleles = parse(Int, parts.nalleles)
                r = readrates(joinpath(root, f))
                out = joinpath(root, replace(f, "rates" => "RNA", ".txt" => ".csv"))
                transitions = get_transitions(G, parts.label)
                write_RNAhistogram(r, transitions, G, R, r[num_rates(transitions, R, S, insertstep)], nalleles, nRNA, outfile=out, ejectnumber=ejectnumber)
            end
        end
    end
end

"""
    tcomponent(model)

return tcomponent of model
"""
tcomponent(model) = typeof(model.components) == TComponents ? model.components : model.components.tcomponents


function make_vector(x, n)
    if !(typeof(x) <: Vector)
        return fill(x, n)
    else
        return x
    end
end


#########
# New functions
#########

function make_observation_dist(d, states, G, R, S, coupling=tuple)
    observations = Vector[]
    if isempty(coupling)
        if eltype(d) <: Distribution
            for s in states
                push!(observations, [d[s] for s in s])
            end
            return states, observations
        else
            units = Vector[]
            for s in eachindex(states)
                push!(observations, [d[s][i] for i in states[s]])
            end
            return states, observations
        end
    else
        if eltype(eltype(d)) <: Distribution
            units = []
            for s in eachindex(states)
                push!(units, [unit_state(i, G, R, S, coupling[1]) for i in states[s]])
                push!(observations, [[d[i] for d in d] for i in states[s]])
            end
            return units, observations
        else
            units = []
            for s in eachindex(states)
                push!(units, [unit_state(i, G, R, S, coupling[1]) for i in states[s]])
                push!(observations, [[d[i] for d in d[s]] for i in states[s]])
            end
            return units, observations
        end
    end

end

function make_trace_datamodel(traces::Vector, interval::Float64, rin, transitions, G, R, S, insertstep, probfn=prob_Gaussian, noiseparams=4, splicetype="", state=true, hierarchical=false, coupling=tuple(), grid=nothing, zeromedian=false)
    data = TraceData{String,String,Tuple}("", "", interval, (traces, [], 0.0, length(traces[1])))
    make_trace_datamodel(data, rin, transitions, G, R, S, insertstep, probfn, noiseparams, splicetype, state, hierarchical, coupling, grid, zeromedian)
end

function make_trace_datamodel(data::TraceData, rin, transitions, G, R, S, insertstep, probfn=prob_Gaussian, noiseparams=4, splicetype="", state=true, hierarchical=false, coupling=tuple(), grid=nothing, zeromedian=false)
    if hierarchical
        h = (2, [8], ())
        method = (Tsit5(), true)
    else
        h = ()
        method = Tsit5()
    end
    if !isempty(coupling)
        model = load_model(data, rin, rin, [1, 2, 3], (), transitions, G, R, S, insertstep, splicetype, 1, 10.0, Int[], 1.0, 0.1, probfn, [ones(Int, noiseparams), ones(Int, noiseparams)], method, h, coupling, grid, zeromedian)
    else
        model = load_model(data, rin, rin, [1, 2, 3], (), transitions, G, R, S, insertstep, splicetype, 1, 10.0, Int[], 1.0, 0.1, probfn, ones(Int, noiseparams), method, h, (), grid, zeromedian)
    end
    return data, model
end


function make_traces_dataframe(ts, td, traces, G::Int, R::Int, S::Int, insertstep::Int, state::Bool, coupling)
    l = maximum(length.(traces))
    data = ["data$i" => [traces[i]; fill(missing, l - length(traces[i]))] for i in eachindex(traces)]
    pred = ["model_mean$i" => [mean.(td[i]); fill(missing, l - length(td[i]))] for i in eachindex(td)]
    predstd = ["model_std$i" => [std.(td[i]); fill(missing, l - length(td[i]))] for i in eachindex(td)]
    cols = [data pred predstd]
    if state
        g, z, zdigits, r = inverse_state(ts, G, R, S, insertstep)
        gs = ["Gstate$i" => [g[i]; fill(missing, l - length(g[i]))] for i in eachindex(g)]
        # tss = ["State$i" => [ts[i]; fill(missing, l - length(g[i]))] for i in eachindex(g)]
        s = ["Rstate$i" => [zdigits[i]; fill(missing, l - length(g[i]))] for i in eachindex(g)]
        # ss = ["Z$i" => [z[i]; fill(missing, l - length(g[i]))] for i in eachindex(g)]
        zs = ["Reporters$i" => [r[i]; fill(missing, l - length(g[i]))] for i in eachindex(g)]
        cols = hcat(cols, [gs s zs])
    end
    # v = state ? [data pred ["state$i" => [mod.(ts[i] .- 1, G) .+ 1; fill(missing, l - length(ts[i]))] for i in eachindex(ts)]] : [data pred]
    # df = DataFrame(["trace$i" => [tp[i]; fill(missing, l - length(tp[i]))] for i in eachindex(tp)])
    DataFrame(permutedims(cols, (2, 1))[:])
end

function make_traces_dataframe(ts, tp, traces, G::Tuple, R::Tuple, S::Tuple, insertstep::Tuple, state::Bool, coupling)
    l = maximum(size.(traces, 1))
    cols = Matrix(undef, length(traces), 0)
    for k in coupling[1]
        data = ["data$i" * "_$k" => [traces[i][:, k]; fill(missing, l - length(traces[i][:, k]))] for i in eachindex(traces)]
        pred = ["model_mean$i" * "_$k" => [[mean(t[k]) for t in tp[i]]; fill(missing, l - length(tp[i]))] for i in eachindex(tp)]
        predstd = ["std_mean$i" * "_$k" => [[std(t[k]) for t in tp[i]]; fill(missing, l - length(tp[i]))] for i in eachindex(tp)]
        cols = hcat(cols, [data pred predstd])
        if state
            index = [[s[k] for s in t] for t in ts]
            g, z, zdigits, r = inverse_state(index, G[k], R[k], S[k], insertstep[k])
            gs = ["Gstate$i" * "_$k" => [g[i]; fill(missing, l - length(g[i]))] for i in eachindex(g)]
            # tss = ["State$i" => [ts[i]; fill(missing, l - length(g[i]))] for i in eachindex(g)]
            s = ["Rstate$i" * "_$k" => [zdigits[i]; fill(missing, l - length(g[i]))] for i in eachindex(g)]
            # ss = ["Z$i" => [z[i]; fill(missing, l - length(g[i]))] for i in eachindex(g)]
            zs = ["Reporters$i" * "_$k" => [r[i]; fill(missing, l - length(g[i]))] for i in eachindex(g)]
            cols = hcat(cols, [gs s zs])
        end
    end
    DataFrame(permutedims(cols, (2, 1))[:])
end

function make_traces_dataframe(traces, interval, rin, transitions, G, R, S, insertstep, probfn=prob_Gaussian, noiseparams=4, splicetype="", state=true, hierarchical=false, coupling=tuple(), grid=nothing, zeromedian=false)
    data = TraceData{String,String,Tuple}("", "", interval, (traces, [], 0.0, length(traces[1])))
    if hierarchical
        h = (2, [8], ())
        method = (Tsit5(), true)
    else
        h = ()
        method = Tsit5()
    end
    if !isempty(coupling)
        model = load_model(data, rin, rin, [1, 2, 3], (), transitions, G, R, S, insertstep, splicetype, 1, 10.0, Int[], 1.0, 0.1, probfn, [ones(Int, noiseparams), ones(Int, noiseparams)], method, h, coupling, grid, zeromedian)
    else
        model = load_model(data, rin, rin, [1, 2, 3], (), transitions, G, R, S, insertstep, splicetype, 1, 10.0, Int[], 1.0, 0.1, probfn, ones(Int, noiseparams), method, h, (), grid, zeromedian)
    end
    ts, d = predict_trace(get_param(model), data, model)
    states, observations = make_observation_dist(d, ts, G, R, S, coupling)
    make_traces_dataframe(states, observations, traces, G, R, S, insertstep, state, coupling)
end

"""
    write_trace_dataframe(outfile, datapath, datacond, interval::Float64, r::Vector, transitions, G, R, S, insertstep, start::Int=1, stop=-1, probfn=prob_Gaussian, noiseparams=4, splicetype=""; state=true, hierarchical=false, coupling=tuple())


"""
function write_trace_dataframe(outfile, datapath, datacond, interval::Float64, r::Vector, transitions, G, R, S, insertstep, start=1, stop=-1, probfn=prob_Gaussian, noiseparams=4, splicetype=""; state=true, hierarchical=false, coupling=tuple(), grid=nothing, zeromedian=false)
    traces = read_tracefiles(datapath, datacond, start, stop)
    traces, maxmedians = zero_median(traces, zeromedian)
    df = make_traces_dataframe(traces, interval, r, transitions, G, R, S, insertstep, probfn, noiseparams, splicetype, state, hierarchical, coupling, grid, zeromedian)
    CSV.write(outfile, df)
end

"""
    write_trace_dataframe(file, datapath::String, datacond, interval, ratetype::String="median", start=1, stop=-1, probfn=prob_Gaussian, noiseparams=4, splicetype=""; hlabel="-h", state=true, coupling=tuple())

TBW
"""
function write_trace_dataframe(file, datapath::String, datacond, interval, ratetype::String="median", start=1, stop=-1, probfn=prob_Gaussian, noiseparams=4, splicetype=""; hlabel="-h", state=true, coupling=tuple(), grid=nothing, zeromedian=false)
    println(file)
    filename = basename(file)
    occursin(hlabel, filename) ? hierarchical = true : hierarchical = false
    parts = fields(filename)
    G, R, S, insertstep = decompose_model(parts.model)
    transitions = get_transitions(G, parts.label)
    r = readrates(file, get_row(ratetype))
    out = replace(file, "rates" => "predictedtraces", ".txt" => ".csv")
    write_trace_dataframe(out, datapath, datacond, interval, r, transitions, G, R, S, insertstep, start, stop, probfn, noiseparams, splicetype, state=state, hierarchical=hierarchical, coupling=coupling, grid=grid, zeromedian=zeromedian)
end


"""
    write_traces(folder, datapath::String, datacond, interval, ratetype::String="median", start=1, stop=-1, probfn=prob_Gaussian, noiseparams=4, splicetype=""; hlabel="-h", state=true, coupling=tuple())

"""
function write_traces(folder, datapath::String, datacond, interval, ratetype::String="median", start=1, stop=-1, probfn=prob_Gaussian, noiseparams=4, splicetype=""; hlabel="-h", state=true, coupling=tuple(), grid=nothing, zeromedian=false)
    for (root, dirs, files) in walkdir(folder)
        for f in files
            # if occursin("rates", f) && occursin(datacond, f) #&& ((!exclude_label && occursin(hlabel, f)) || exclude_label && !occursin(hlabel, f))
            if occursin("rates", f) && (occursin("tracejoint", f) || (isempty(coupling) && occursin(datacond, f)))
                write_trace_dataframe(joinpath(root, f), datapath, datacond, interval, ratetype, start, stop, probfn, noiseparams, splicetype, hlabel=hlabel, state=state, coupling=coupling, grid=grid, zeromedian=zeromedian)
            end
        end
    end
end
"""
    write_traces_folder(folder, datafolder::Vector, datacond, interval, ratetype::String="median", start=1, stop=-1, probfn=prob_Gaussian, noiseparams=4, splicetype=""; state=true, coupling=tuple())

TBW
"""
function write_traces_coupling(folder, datapath, datacond, interval, G=(3, 3), R=(3, 3), sources=1:3, targets=1:5, ratetype::String="median", start=1, stop=-1, probfn=prob_Gaussian, noiseparams=4, splicetype=""; hlabel="-h", state=true, pattern="gene", zeromedian=true)
    for (root, dirs, files) in walkdir(folder)
        for f in files
            for target in targets
                for source in sources
                    # if occursin("rates", f) && occursin(datacond, f) #&& ((!exclude_label && occursin(hlabel, f)) || exclude_label && !occursin(hlabel, f))
                    if occursin("rates", f) && occursin("$pattern$source$target", f)
                        write_trace_dataframe(joinpath(root, f), datapath, datacond, interval, ratetype, start, stop, probfn, noiseparams, splicetype, hlabel=hlabel, state=state, coupling=((1, 2), (tuple(), tuple(1)), (source, 0), (0, target), 1), zeromedian=zeromedian)
                    end
                end
                if occursin("rates", f) && occursin("R$target", f)
                    coupling = ((1, 2), (tuple(), tuple(1)), (collect(G[1]+1:G[1]+R[1]), 0), (0, target), 1)
                    write_trace_dataframe(joinpath(root, f), datapath, datacond, interval, ratetype, start, stop, probfn, noiseparams, splicetype, hlabel=hlabel, state=state, coupling=coupling, zeromedian=zeromedian)
                end
            end
        end
    end
end

########### Parallelized functions ###########
# """
#     write_traces_coupling_parallel(folder, datapath, datacond, interval, G=(3, 3), R=(3, 3), 
#                                    sources=1:3, targets=1:5, ratetype::String="median", 
#                                    start=1, stop=-1, probfn=prob_Gaussian, noiseparams=4, 
#                                    splicetype=""; hlabel="-h", state=true, pattern="gene")

# Parallelized version of write_traces_coupling using Julia's multi-threading.
# """
# function write_traces_coupling_parallel(folder, datapath, datacond, interval, G=(3, 3), R=(3, 3),
#     sources=1:3, targets=1:5, ratetype::String="median",
#     start=1, stop=-1, probfn=prob_Gaussian, noiseparams=4,
#     splicetype=""; hlabel="-h", state=true, pattern="gene")

#     # Collect all tasks that need to be executed
#     tasks = []

#     for (root, dirs, files) in walkdir(folder)
#         for f in files
#             if occursin("rates", f)
#                 # Tasks for pattern$source$target files
#                 for target in targets
#                     for source in sources
#                         if occursin("$pattern$source$target", f)
#                             coupling = ((1, 2), (tuple(), tuple(1)), (source, 0), (0, target), 1)
#                             task_args = (
#                                 joinpath(root, f), datapath, datacond, interval, ratetype,
#                                 start, stop, probfn, noiseparams, splicetype,
#                                 hlabel, state, coupling
#                             )
#                             push!(tasks, task_args)
#                         end
#                     end

#                     # Tasks for R$target files
#                     if occursin("R$target", f)
#                         coupling = ((1, 2), (tuple(), tuple(1)), (collect(G[1]+1:G[1]+R[1]), 0), (0, target), 1)
#                         task_args = (
#                             joinpath(root, f), datapath, datacond, interval, ratetype,
#                             start, stop, probfn, noiseparams, splicetype,
#                             hlabel, state, coupling
#                         )
#                         push!(tasks, task_args)
#                     end
#                 end
#             end
#         end
#     end

#     # Process each task in parallel using threads
#     @info "Processing $(length(tasks)) files in parallel using $(Threads.nthreads()) threads"

#     Threads.@threads for task_args in tasks
#         file_path, datapath, datacond, interval, ratetype, start, stop, probfn,
#         noiseparams, splicetype, hlabel, state, coupling = task_args

#         try
#             write_trace_dataframe(file_path, datapath, datacond, interval, ratetype,
#                 start, stop, probfn, noiseparams, splicetype,
#                 hlabel=hlabel, state=state, coupling=coupling)
#             @info "Successfully processed: $(basename(file_path))"
#         catch e
#             @error "Error processing file: $(basename(file_path))" exception = (e, catch_backtrace())
#         end
#     end
# end

# Alternative implementation using Threads.@spawn for more dynamic scheduling
"""
    write_traces_coupling_spawn(folder, datapath, datacond, interval, G=(3, 3), R=(3, 3),
                               sources=1:3, targets=1:5, ratetype::String="median",
                               start=1, stop=-1, probfn=prob_Gaussian, noiseparams=4,
                               splicetype=""; hlabel="-h", state=true, pattern="gene")

Parallelized version of write_traces_coupling using Julia's task-based parallelism with @spawn.
This provides more dynamic scheduling which can be beneficial for workloads with varying execution times.
"""
function write_traces_coupling_spawn(folder, datapath, datacond, interval, G=(3, 3), R=(3, 3),
    sources=1:3, targets=1:5, ratetype::String="median",
    start=1, stop=-1, probfn=prob_Gaussian, noiseparams=4,
    splicetype=""; hlabel="-h", state=true, pattern="gene")

    # Collect all tasks that need to be executed
    tasks = []

    for (root, dirs, files) in walkdir(folder)
        for f in files
            if occursin("rates", f)
                # Tasks for pattern$source$target files
                for target in targets
                    for source in sources
                        if occursin("$pattern$source$target", f)
                            file_path = joinpath(root, f)
                            coupling = ((1, 2), (tuple(), tuple(1)), (source, 0), (0, target), 1)

                            # Create a task for this file
                            t = @spawn begin
                                try
                                    write_trace_dataframe(file_path, datapath, datacond, interval, ratetype,
                                        start, stop, probfn, noiseparams, splicetype,
                                        hlabel=hlabel, state=state, coupling=coupling)
                                    @info "Successfully processed: $(basename(file_path))"
                                catch e
                                    @error "Error processing file: $(basename(file_path))" exception = (e, catch_backtrace())
                                end
                            end

                            push!(tasks, t)
                        end
                    end

                    # Tasks for R$target files
                    if occursin("R$target", f)
                        file_path = joinpath(root, f)
                        coupling = ((1, 2), (tuple(), tuple(1)), (collect(G[1]+1:G[1]+R[1]), 0), (0, target), 1)

                        # Create a task for this file
                        t = @spawn begin
                            try
                                write_trace_dataframe(file_path, datapath, datacond, interval, ratetype,
                                    start, stop, probfn, noiseparams, splicetype,
                                    hlabel=hlabel, state=state, coupling=coupling)
                                @info "Successfully processed: $(basename(file_path))"
                            catch e
                                @error "Error processing file: $(basename(file_path))" exception = (e, catch_backtrace())
                            end
                        end

                        push!(tasks, t)
                    end
                end
            end
        end
    end

    @info "Scheduled $(length(tasks)) files for processing using $(Threads.nthreads()) threads"

    # Wait for all tasks to complete
    for t in tasks
        fetch(t)
    end

    @info "All files processed"
end

# Usage example:
# First set the number of threads via environment variable or command line:
# $ julia -t 8 your_script.jl
# Or within Julia:
# using Base.Threads
# @show nthreads()
#
# Then call either of the parallel functions:
# write_traces_coupling_parallel(folder, datapath, datacond, interval)
# write_traces_coupling_spawn(folder, datapath, datacond, interval)



function extract_source_target(pattern::String, filepath::String)
    # Get filename from path
    filename = basename(filepath)

    # Escape special characters in pattern and build regex
    escaped_pattern = replace(pattern, r"([.*+?^\$()[]{}|\\])" => s"\\\1")
    regex = Regex("$(escaped_pattern)(\\d|R)(\\d)")

    # Find pattern## in filename
    m = match(regex, filename)
    if m !== nothing
        source = m[1] == "R" ? "R" : parse(Int, m[1])
        target = parse(Int, m[2])
        return source, target
    end
    return nothing
end

function write_cov_file(file, transitions=(([1, 2], [2, 1], [2, 3], [3, 2]), ([1, 2], [2, 1], [2, 3], [3, 2])), G=(3, 3), R=(3, 3), S=(1, 0), insertstep=(1, 1), interval=1.0, pattern="gene", lags=collect(0:10:1000), probfn=prob_Gaussian, ratetype="ml")
    println(file)
    r = readrates(file, get_row(ratetype))
    source, target = extract_source_target(pattern, file)
    (source == "R") && (source = collect(G[1]+1:G[1]+R[1]))
    coupling = ((1, 2), (tuple(), tuple(1)), (source, 0), (0, target), 1)
    covariance_functions(r, transitions, G, R, S, insertstep, interval, probfn, coupling, lags)
end

function write_cov(folder, transitions=(([1, 2], [2, 1], [2, 3], [3, 2]), ([1, 2], [2, 1], [2, 3], [3, 2])), G=(3, 3), R=(3, 3), S=(0, 0), insertstep=(1, 1), interval=1.0, pattern="gene", lags=collect(0:10:1000), probfn=prob_Gaussian, ratetype="ml")
    for (root, dirs, files) in walkdir(folder)
        for f in files
            if occursin("rates", f) && occursin("tracejoint", f)
                file = joinpath(root, f)
                ac1, ac2, cc, ac1f, ac2f, ccf, tau, _, _, _, _ = write_cov_file(file, transitions, G, R, S, insertstep, interval, pattern, lags, probfn, ratetype)
                out = replace(file, "rates" => "crosscovariance", ".txt" => ".csv")
                CSV.write(out, DataFrame(tau=tau, crosscovariance=cc, autocovariance1=ac1, autocovariance2=ac2, cc_normalized=ccf, ac1_normalized=ac1f, ac2_normalized=ac2f))
            end
        end
    end
end

function write_residency_G_allgenes(fileout::String, filein::String, G, header)
    rates = read_all_rates_csv(filein, header)
    n = G - 1
    f = open(fileout, "w")
    writedlm(f, ["gene" collect(0:n)'], ',')
    for r in eachrow(rates)
        writedlm(f, [r[1] residenceprob_G(r[2:2*n+1], n + 1)], ',')
    end
    close(f)
end

function write_residency_G(file, ratetype="median", transitions=(([1, 2], [2, 1], [2, 3], [3, 2]), ([1, 2], [2, 1], [2, 3], [3, 2])), nnoiseparams=4)
    println(file)
    r = readrates(file, get_row(ratetype))
    parts = fields(basename(file))
    G, R, S, insertstep = decompose_model(parts.model)
    out = replace(file, "rates" => "residencyG", ".txt" => ".csv")
    if G isa Int
        CSV.write(out, DataFrame(residenceprob_G_dataframe(r, G)))
    else
        if G[1] == 2 && G[2] == 2
            transitions = (([1, 2], [2, 1]), ([1, 2], [2, 1]))
        end
        nrates = num_rates(transitions, R, S, insertstep) .+ nnoiseparams
        CSV.write(out, DataFrame(residenceprob_G_dataframe(r, G, nrates)))
    end
end

function write_residency_G_folder(folder)
    for (root, dirs, files) in walkdir(folder)
        for f in files
            if occursin("rates", f)
                file = joinpath(root, f)
                write_residency_G(file)
            end
        end
    end

end

"""
    make_trace_histogram(datapath, datacond, start=1, stop=-1)

TBW
"""
function make_trace_histogram(datapath, datacond, start=1, stop=-1)
    traces = read_tracefiles(datapath, datacond, start, stop)
    ft = reduce(vcat, traces)
    h = histogram(ft, normalize=true)
    return ft, h
end





"""
    simulate_trials(r, transitions, G, R, S, insertstep, ntrials, trial_time=720.)

Simulate multiple trials of a gene expression system and compare with theoretical predictions.

# Arguments
- `r`: Rate parameters vector
- `transitions`: State transition matrix
- `G`: Gene states matrix
- `R`: RNA states matrix
- `S`: Initial state vector
- `insertstep`: Time step for trajectory recording
- `ntrials`: Number of simulation trials
- `trial_time`: Total simulation time (default: 720.0)

# Returns
- Tuple containing:
    1. `(linf_norm, l2_norm)`: L-infinity and L2 norms between theory and simulation
    2. `cc_theory`: Theoretical steady-state value
    3. `cc_simulated`: Array of maximum values from each trial
    4. `running_means`: Running mean of simulated values
    5. `running_errors`: Running standard errors of the mean

# Example
```julia
ntrials = 100
(norms, theory, sims, means, errors) = simulate_trials(r, trans, G, R, S, 0.1, ntrials)
plot(1:ntrials, errors, ylabel="Standard Error", xlabel="Number of Trials")
"""

function simulate_trials(r::Vector, transitions::Tuple, G, R, S, insertstep, coupling, ntrials, trial_time=720.0, lag=60, stride=1, probfn=prob_Gaussian)
    _, _, _, _, ac1, ac2, cc_theory, _ = StochasticGene.covariance_functions(r, transitions, G, R, S, insertstep, 1.0, probfn, coupling, collect(0:stride:lag))
    simulate_trials(ac1, ac2, cc_theory, r, transitions, G, R, S, insertstep, coupling, ntrials, trial_time, lag, stride)
end

function simulate_trials(ac1_theory, ac2_theory, cc_theory::Vector, r::Vector, transitions::Tuple, G, R, S, insertstep, coupling, ntrials, trial_time::Float64=720.0, lag=60, stride=1, col=2)
    lags = collect(-lag:stride:lag)
    cc_mean = zeros(length(lags))
    cc_var = zeros(length(lags))
    lags_ac = collect(0:stride:lag)
    ac1_mean = zeros(length(lags_ac))
    ac1_var = zeros(length(lags_ac))
    ac2_mean = zeros(length(lags_ac))
    ac2_var = zeros(length(lags_ac))
    linf_norm = Float64[]
    l2_norm = Float64[]
    ac1tuple = (0, ac1_mean, ac1_var)
    ac2tuple = (0, ac2_mean, ac2_var)
    vartuple = (0, cc_mean, cc_var)
    for n in 1:ntrials
        t = simulate_trace_vector(r, transitions, G, R, S, insertstep, coupling, 1.0, trial_time, 1, col=col)
        vartuple = var_update(vartuple, StatsBase.crosscov(t[1][:, 1], t[1][:, 2], lags, demean=true))
        ac1tuple = var_update(ac1tuple, StatsBase.autocov(t[1][:, 1], collect(0:stride:lag), demean=true))
        ac2tuple = var_update(ac2tuple, StatsBase.autocov(t[1][:, 2], collect(0:stride:lag), demean=true))
        push!(linf_norm, maximum(abs.(cc_theory .- vartuple[2])))
        push!(l2_norm, sqrt(sum((cc_theory .- vartuple[2]) .^ 2)))
    end
    counts, cc_mean, cc_var = vartuple
    counts1, ac1_mean, ac1_var = ac1tuple
    counts2, ac2_mean, ac2_var = ac2tuple
    return linf_norm, l2_norm, cc_theory, cc_mean, sqrt.(cc_var ./ (counts - 1) ./ counts), lags, ac1_mean, ac1_theory, ac2_mean, ac2_theory
end

function data_covariance(traces, lags)
    ac1 = zeros(length(lags))
    ac2 = zeros(length(lags))
    cov = zeros(length(lags))
    for t in traces
        cov += StatsBase.crosscov(t[:, 1], t[:, 2], lags, demean=true)
        ac1 += StatsBase.crosscov(t[:, 1], t[:, 1], lags, demean=true)
        ac2 += StatsBase.crosscov(t[:, 2], t[:, 2], lags, demean=true)
    end
    ac1 / length(traces), ac2 / length(traces), cov / length(traces), lags
end

function reset_noise!(r, transitions, G::Tuple, R, S, insertstep, num_noiseparams)
    n = 0
    for i in eachindex(G)
        n += num_rates(transitions[i], R[i], S[i], insertstep[i])
        r[n+1:n+num_noiseparams] .= [20, 0.1, 80, 0.1]
        n += num_noiseparams
    end
end

function delete_noise!(r, transitions, G::Tuple, R, S, insertstep, num_noiseparams)
    n = 0
    for i in eachindex(G)
        n += num_rates(transitions[i], R[i], S[i], insertstep[i])
        deleteat!(r, n+1:n+num_noiseparams)
        # n += num_noiseparams
    end
end

"""
    plot_traces(df::DataFrame; 
               cols::Vector{Int}=collect(1:3),  # Convert UnitRange to Vector
               with_errors::Bool=true,
               normalize::Bool=true,
               title::String="Traces",
               alpha::Float64=0.3)

Plot multiple traces in separate panels for data and model.

# Arguments
- `df::DataFrame`: DataFrame containing data/model columns
- `cols::Vector{Int}`: Column numbers to plot (default: [1,2,3])
- `with_errors::Bool`: Include error ribbons (default: true)
- `normalize::Bool`: Normalize data (default: true)
- `title::String`: Plot title (default: "Traces")
- `alpha::Float64`: Transparency for ribbons (default: 0.3)

# Returns
- `Plots.Plot`: Combined plot object
"""
function plot_traces(df::DataFrame;
    cols::Vector{Int}=collect(1:2),  # Changed this line
    with_errors::Bool=false,
    normalize::Bool=false,
    title::String="Traces",
    alpha::Float64=0.3)

    n = length(cols)
    p = plot(layout=(n, 1), legend=:topleft, title=title)

    for (j, i) in enumerate(cols)
        data_col = Symbol("data$i")
        mean_col = Symbol("model_mean$i")
        std_col = Symbol("model_std$i")

        # Get data, filtering out missing values
        y_data = df[!, data_col]
        y_model = df[!, mean_col]
        y_std = hasproperty(df, std_col) ? df[!, std_col] : nothing

        # Filter out missing values
        valid_indices = .!ismissing.(y_data) .& .!ismissing.(y_model)
        y_data = y_data[valid_indices]
        y_model = y_model[valid_indices]
        y_std = y_std !== nothing ? y_std[valid_indices] : nothing

        # Normalize if needed
        if normalize
            y_data ./= maximum(y_data)
            y_model ./= maximum(y_model)
            y_std = y_std !== nothing ? y_std ./ maximum(y_model) : nothing
        end

        # Plot data points in the first row
        scatter!(p[j, 1], y_data, label="Data $i", markersize=3)

        # Plot model with optional error ribbon in the second row
        if with_errors && y_std !== nothing
            plot!(p[j, 1], y_model, err=y_std, label="Model $i", alpha=alpha)
        else
            plot!(p[j, 1], y_model, label="Model $i")
        end
    end

    # xlabel!(p[1, :], "Time point")
    # ylabel!(p[1, :], normalize ? "Normalized value" : "Value")
    # xlabel!(p[2, :], "Time point")
    # ylabel!(p[2, :], normalize ? "Normalized value" : "Value")

    return p
end


"""
    plot_histogram(ratefile::String, datapath; root=".", row=2)

    plot_histogram()
    plot_histogram(ratefile::String,datapath;fish=false,root=".",row=2)
    
    functions to plot data and model predicted histograms
    
"""
function plot_histogram(ratefile::String, datapath; root=".", row=2)
    fish = false
    r = readrow(ratefile, row)
    println("r: ", r)
    parts = fields(ratefile)
    label = parts.label
    cond = parts.cond
    G = parts.model
    data = data_rna(parts.gene, parts.cond, datapath, fish, parts.label, root)
    model = model_rna(r, r, parse(Int, parts.model), parse(Int, parts.nalleles), 0.01, [], (), 0)
    # model_rna(r,[rateprior[i]],G,nalleles,cv,[fittedparam[i]],fixedeffects,0)
    plot_histogram(data, model)
    return data, model
end

function plot_histogram(gene::String, cell::String, G::Int, cond::String, ratefile::String, datapath::String, root::String=".")
    fish = false
    rates = readdlm(ratefile, ',', header=true)
    r = rates[1][findfirst(rates[1][:, 1] .== gene)[1], 2:end]
    data = data_rna(gene, cond, datapath, fish, "label", root)
    nalleles = alleles(gene, cell, root)
    model = model_rna(r, [], G, nalleles, 0.01, [], (), 0)
    println("typeof(model): ", typeof(model))
    println("typeof(data): ", typeof(data))
    m = plot_histogram(data, model)
    return m, data, model
end

function plot_histogram(gene::String, cell::String, G::String, cond::String, label::String, ratefolder::String, datapath::String, nsets::Int, root::String, fittedparam=[1], verbose=false)
    fish = false
    data = data_rna(gene, cond, datapath, fish, label, root)
    model = model_rna(gene, cell, G, fish, 0.01, fittedparam, (), label, ratefolder, nsets, root, data, verbose)
    m = plot_histogram(data, model)
    return m, data, model
end

function plot_histogram(data::AbstractRNAData{Array{Array,1}}, model)
    h = predictedarray(model.rates, data, model)
    figure(data.gene)
    for i in eachindex(h)
        plot(h[i])
        plot(normalize_histogram(data.histRNA[i]))
        savefig(string(i))
    end
    return h
end
function plot_histogram(data::AbstractRNAData{Array{Float64,1}}, model)
    h = predictedfn(get_param(model), data, model)
    plt = plot(h)
    plot!(plt, normalize_histogram(data.histRNA))
    display(plt)
    return h
end

function plot_histogram(data::RNAOnOffData, model::AbstractGeneTransitionModel, filename="")
    h = predictedarray(model.rates, data, model)
    plt1 = plot(data.bins, h[1])
    plot!(plt1, data.bins, normalize_histogram(data.OFF))
    plt2 = plot(data.bins, h[2])
    plot!(plt2, data.bins, normalize_histogram(data.ON))
    plt3 = plot(h[3])
    plot!(plt3, normalize_histogram(data.histRNA))
    plt = plot(plt1, plt2, plt3, layout=(3, 1))
    display(plt)
    if ~isempty(filename)
        savefig(filename)
    end
    return h
end

function plot_histogram(data::TraceRNAData, model::AbstractGeneTransitionModel, filename="")
    M = make_mat_M(model.components.mcomponents, model.rates)
    h = steady_state(M, model.components.mcomponents.nT, model.nalleles, data.nRNA)
    plt = plot(h)
    plot!(plt, normalize_histogram(data.histRNA))
    display(plt)
    if ~isempty(filename)
        savefig(filename)
    end
    return h
end

function plot_histogram(data::RNAData{T1,T2}, model::AbstractGMmodel, save=false) where {T1<:Array,T2<:Array}
    m = predictedarray(model.rates, data, model)
    for i in eachindex(m)
        plt = plot(m[i])
        plot!(normalize_histogram(data.histRNA[i]), show=true)
        if save
            savefig()
        else
            display(plt)
        end
    end
    println("deviance: ", deviance(data, model))
    println("loglikelihood: ", loglikelihood(get_param(model), data, model)[1])
    return m
end

function plot_histogram(data::RNAData, model::AbstractGMmodel)
    h = predictedfn(get_param(model), data, model)
    plt = plot(h)
    plot!(normalize_histogram(data.histRNA))
    display(plt)
    return h
end

# Utility function for splitting conditions
split_conditions(cond::AbstractString, multicond::Bool) = multicond ? split(cond, "-") : [cond]

