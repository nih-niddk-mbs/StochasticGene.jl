# This file is part of StochasticGene.jl   

# analysis.jl

"""
    make_dataframes(resultfolder::String, datapath::String, assemble=true, multicond=false, datatype="rna")

Create and assemble dataframes from model fitting results.

# Arguments
- `resultfolder::String`: Path to folder containing result files
- `datapath::String`: Path to folder containing input data files
- `assemble::Bool=true`: Whether to assemble results into summary files
- `multicond::Bool=false`: Whether to handle multiple conditions
- `datatype::String="rna"`: Type of data ("rna", "rnacount", etc.)

# Returns
- `Vector{Vector}`: Nested vector of tuples containing (filename, DataFrame) pairs

# Notes
- Automatically assembles results if assemble=true
- Processes rate summary files from the result folder
- Groups files by label and model type
- Creates summary files for each label-model combination
- Handles both single and multiple condition datasets
- Supports different data types (RNA, RNA count, etc.)

# Examples
```julia
# Basic usage
dfs = make_dataframes("results/", "data/")

# With multiple conditions
dfs = make_dataframes("results/", "data/", assemble=true, multicond=true)

# For RNA count data
dfs = make_dataframes("results/", "data/", datatype="rnacount")
```
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

Convert a rate file path to its corresponding statistics file path.

# Arguments
- `ratefile::String`: Path to a rate file

# Returns
- `String`: Path to the corresponding statistics file

# Notes
- Replaces "rates_" with "stats_" in the filename
- Assumes rate files and stat files follow the same naming convention
- Used to find corresponding statistics files for rate files

# Examples
```julia
ratefile = "rates_gene_condition_3401_2.txt"
statfile = statfile_from_ratefile(ratefile)
# Returns: "stats_gene_condition_3401_2.txt"
```
"""
statfile_from_ratefile(ratefile) = replace(ratefile, "rates_" => "stats_")

"""
    make_dataframe(ratefile::String, datapath::String, multicond::Bool, datatype::String)

Create a DataFrame from a rate file with additional data and statistics.

# Arguments
- `ratefile::String`: Path to the rate file
- `datapath::String`: Path to the data directory
- `multicond::Bool`: Whether to handle multiple conditions
- `datatype::String`: Type of data ("rna", "rnacount", etc.)

# Returns
- `DataFrame`: DataFrame containing rate data, statistics, and moments

# Notes
- Reads rate data and corresponding statistics
- Joins rate and statistics dataframes
- Adds model information based on filename parsing
- Handles multiple conditions if specified
- Adds expression moments (mean, variance, third moment)
- Only processes single-state models (G models)
- Skips multi-state models with a warning message

# Examples
```julia
# Basic usage
df = make_dataframe("rates_gene_condition_3401_2.txt", "data/", false, "rna")

# With multiple conditions
df = make_dataframe("rates_gene_condition_3401_2.txt", "data/", true, "rnacount")
```
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

Parse a model name string to extract model parameters.

# Arguments
- `name::String`: Model name string (e.g., "3401", "2121")

# Returns
- `Vector{Int}`: Vector of model parameters (G, R, S, insertstep)

# Notes
- Converts string to integer and extracts digits
- For single-digit models: returns single integer
- For multi-digit models: returns vector of digits
- Assumes model parameters are encoded as consecutive digits
- Used to extract G, R, S, and insertstep values from model names

# Examples
```julia
# Single state model
parse_model("3401")  # Returns: [3, 4, 0, 1]

# Multi-state model
parse_model("2121")  # Returns: [2, 1, 2, 1]
```
"""
function parse_model(name::String)
    d = parse(Int, name)
    d > 9 && (d = digits(d))
    d
end

"""
    augment_dataframe(df, resultfolder)

Augment a DataFrame with additional model-specific information and statistics.

# Arguments
- `df::DataFrame`: Input DataFrame containing model results
- `resultfolder::String`: Path to folder containing result files

# Returns
- `DataFrame`: Augmented DataFrame with additional columns

# Notes
- Adds burst size information for models with G > 1
- Adds model moments for G = 2 models
- Adds fit measures (AIC, BIC, etc.) for the model
- Adds residence probabilities for gene states
- Creates a copy of the input DataFrame to avoid modifying original
- Automatically determines which augmentations to apply based on model type

# Examples
```julia
# Augment a DataFrame with all available information
df_augmented = augment_dataframe(df, "results/")

# The augmented DataFrame will contain additional columns like:
# - BurstMean, BurstSD, BurstMedian, BurstMAD (for G > 1)
# - Model_Expression, Model_Variance (for G = 2)
# - AIC, BIC, WAIC, etc. (fit measures)
# - ProbG0, ProbG1, etc. (residence probabilities)
```
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

Create a DataFrame from measure summary files for a specific model.

# Arguments
- `resultfolder::String`: Path to folder containing result files
- `G::String`: Model identifier (e.g., "2", "3")

# Returns
- `DataFrame`: Combined DataFrame from all measure files for the specified model

# Notes
- Searches for measure summary files in the result folder
- Filters files by the specified model identifier
- Adds condition information if available
- Combines all matching files into a single DataFrame
- Uses stack_dataframe to merge multiple files

# Examples
```julia
# Create measure DataFrame for G=2 model
df_measures = make_measure_df("results/", "2")

# The resulting DataFrame contains fit measures like:
# - AIC, BIC, WAIC
# - Deviance, LogLikelihood
# - Gene and Condition columns
```
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

Determine the columns to use for joining DataFrames based on condition availability.

# Arguments
- `d::DataFrame`: Input DataFrame

# Returns
- `Vector{Symbol}`: Vector of column symbols to use for joining

# Notes
- Returns [:Gene] if Condition column is missing or has missing values
- Returns [:Gene, :Condition] if Condition column is available
- Used to determine appropriate join keys for merging DataFrames
- Ensures consistent joining behavior across different data structures

# Examples
```julia
# DataFrame with conditions
df_with_cond = DataFrame(Gene=["gene1", "gene2"], Condition=["ctrl", "ctrl"])
join_cols(df_with_cond)  # Returns: [:Gene, :Condition]

# DataFrame without conditions
df_no_cond = DataFrame(Gene=["gene1", "gene2"])
join_cols(df_no_cond)    # Returns: [:Gene]
```
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

Add fit measures to a DataFrame by joining with measure summary files.

# Arguments
- `df::DataFrame`: Input DataFrame containing model results
- `resultfolder::String`: Path to folder containing result files
- `G`: Model identifier (can be String or other type)

# Returns
- `DataFrame`: DataFrame with added measure columns

# Notes
- Creates measure DataFrame for the specified model
- Joins measure data with input DataFrame using appropriate columns
- Automatically determines join columns based on condition availability
- Adds fit measures like AIC, BIC, WAIC, Deviance, LogLikelihood
- Preserves all original columns from input DataFrame

# Examples
```julia
# Add fit measures to a DataFrame
df_with_measures = add_measures(df, "results/", "2")

# The resulting DataFrame will have additional columns:
# - AIC, BIC, WAIC
# - Deviance, LogLikelihood
# - Other model fit statistics
```
"""
function add_measures(df, resultfolder::String, G)
    dm = make_measure_df(resultfolder, G)
    leftjoin(df, dm, on=join_cols(df))
end

"""
    add_mean!(df::DataFrame, datapath)

Add mean expression values to a DataFrame by calculating from RNA histograms.

# Arguments
- `df::DataFrame`: Input DataFrame (modified in place)
- `datapath::String`: Path to data directory containing RNA histograms

# Returns
- `Nothing`: Modifies the DataFrame in place

# Notes
- Calculates mean expression for each gene-condition combination
- Uses RNA histogram data from the specified data path
- Adds an :Expression column to the DataFrame
- Assumes gene names are in :Gene column and conditions in :Condition column
- Uses mean_histogram function to calculate means from histograms

# Examples
```julia
# Add mean expression values
add_mean!(df, "data/")

# The DataFrame now contains an :Expression column with mean values
# calculated from RNA histogram data
```
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
    add_moments!(df::DataFrame, datapath, datatype::String)

Add expression moments (mean, variance, third moment) to a DataFrame.

# Arguments
- `df::DataFrame`: Input DataFrame (modified in place)
- `datapath::String`: Path to data directory
- `datatype::String`: Type of data ("rna" or "rnacount")

# Returns
- `Nothing`: Modifies the DataFrame in place

# Notes
- Adds :Expression, :Variance, and :ThirdMoment columns
- For "rna" datatype: uses RNA histogram data
- For "rnacount" datatype: uses RNA count data
- Calculates moments for each gene-condition combination
- Uses appropriate moment calculation functions based on data type
- Assumes gene names are in :Gene column and conditions in :Condition column

# Examples
```julia
# Add moments for RNA histogram data
add_moments!(df, "data/", "rna")

# Add moments for RNA count data
add_moments!(df, "data/", "rnacount")

# The DataFrame now contains:
# - :Expression (mean)
# - :Variance (variance)
# - :ThirdMoment (third moment)
```
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

Add model-predicted moments (mean and variance) to a DataFrame for G=2 models.

# Arguments
- `df::DataFrame`: Input DataFrame (modified in place)

# Returns
- `Nothing`: Modifies the DataFrame in place

# Notes
- Adds :Model_Expression and :Model_Variance columns
- Only applies to G=2 models (two-state gene models)
- Uses model2_mean and model2_variance functions
- Requires columns: Rate12, Rate21, Eject, Decay, Nalleles
- Calculates theoretical moments based on fitted model parameters
- Useful for comparing data moments with model predictions

# Examples
```julia
# Add model moments for G=2 model results
add_modelmoments!(df)

# The DataFrame now contains:
# - :Model_Expression (theoretical mean)
# - :Model_Variance (theoretical variance)

# These can be compared with :Expression and :Variance from data
```
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

Add burst size information to a DataFrame by joining with burst data.

# Arguments
- `df::DataFrame`: Input DataFrame to augment
- `db::DataFrame`: DataFrame containing burst size data
- `cols::Vector{Symbol}`: Column names to use for joining

# Returns
- `DataFrame`: DataFrame with added burst size columns

# Notes
- Joins burst data based on specified columns
- Adds columns: BurstMean, BurstSD, BurstMedian, BurstMAD
- Handles both single and multi-condition datasets
- Uses left join to preserve all rows from input DataFrame
- Automatically includes Gene and Condition columns if available

# Examples
```julia
# Add burst size data
df_with_burst = add_burstsize(df, burst_data, [:Gene, :Condition])

# The resulting DataFrame will have additional columns:
# - BurstMean, BurstSD, BurstMedian, BurstMAD
```
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

Create a DataFrame from burst summary files in a result folder.

# Arguments
- `resultfolder::String`: Path to folder containing result files

# Returns
- `DataFrame`: Combined DataFrame from all burst summary files

# Notes
- Searches for burst summary files in the result folder
- Parses filename components to extract condition information
- Adds condition information if available in filename
- Combines all burst files into a single DataFrame
- Uses stack_dataframe to merge multiple files

# Examples
```julia
# Create burst DataFrame from result folder
burst_df = make_burst_df("results/")

# The resulting DataFrame contains:
# - Gene, Condition columns
# - BurstMean, BurstSD, BurstMedian, BurstMAD columns
```
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

Add a timestamp column to a DataFrame.

# Arguments
- `df::DataFrame`: Input DataFrame (modified in place)
- `timestamp`: Timestamp value to add (can be any type)

# Returns
- `Nothing`: Modifies the DataFrame in place

# Notes
- Adds a :Time column to the DataFrame
- Fills all rows with the specified timestamp value
- Useful for tracking when data was processed or analyzed
- Modifies the DataFrame in place for efficiency

# Examples
```julia
# Add current timestamp
add_time!(df, now())

# Add custom timestamp
add_time!(df, "2024-01-15")

# Add numeric timestamp
add_time!(df, 1234567890)
```
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
    separate_dataframe(df, G, cond, multicond)

Separate a DataFrame into multiple DataFrames based on conditions.

# Arguments
- `df::DataFrame`: Input DataFrame
- `G::Int`: Number of gene states
- `cond::String`: Condition string
- `multicond::Bool`: Whether to handle multiple conditions

# Returns
- `Vector{DataFrame}`: Vector of separated DataFrames

# Notes
- Splits condition string if multicond=true
- Creates separate DataFrame for each condition
- Extracts appropriate columns for each condition based on G
- Adds Condition column to each DataFrame
- Renames columns if multicond=true by removing suffixes
- Assumes column structure follows standard format

# Examples
```julia
# Separate DataFrame with multiple conditions
dfs = separate_dataframe(df, 2, "control-treatment", true)

# Returns vector of DataFrames, one for each condition
```
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

"""
    add_Zscore_class!(df, threshold=2.0)

Add Z-score classification columns to a DataFrame based on a threshold.

# Arguments
- `df::DataFrame`: Input DataFrame (modified in place)
- `threshold::Float64=2.0`: Threshold for classifying Z-scores

# Returns
- `Nothing`: Modifies the DataFrame in place

# Notes
- Adds classification columns: :On, :Off, :Eject, :OnOffEject
- Classifies Z-scores as "Up", "Down", or "Null" based on threshold
- Requires existing Z-score columns: :ZRate01, :ZRate10, :ZEject
- :OnOffEject combines all three classifications
- Uses classify_Zscore function for individual classifications
- Modifies the DataFrame in place for efficiency

# Examples
```julia
# Add Z-score classifications with default threshold
add_Zscore_class!(df)

# Add Z-score classifications with custom threshold
add_Zscore_class!(df, 1.5)

# The DataFrame now contains:
# - :On, :Off, :Eject (individual classifications)
# - :OnOffEject (combined classification)
```
"""
function add_Zscore_class!(df, threshold=2.0)
    insertcols!(df, :On => classify_Zscore.(df.ZRate01, threshold))
    insertcols!(df, :Off => classify_Zscore.(df.ZRate10, threshold))
    insertcols!(df, :Eject => classify_Zscore.(df.ZEject, threshold))
    insertcols!(df, :OnOffEject => df.On .* df.Off .* df.Eject)
end

"""
    classify_Zscore(Z, threshold)

Classify a Z-score value based on a threshold.

# Arguments
- `Z::Float64`: Z-score value to classify
- `threshold::Float64`: Threshold for classification

# Returns
- `String`: Classification ("Up", "Down", or "Null")

# Notes
- Returns "Up" if Z > threshold
- Returns "Down" if Z < -threshold
- Returns "Null" if -threshold ≤ Z ≤ threshold
- Used for categorizing changes in gene expression parameters
- Typically used with Z-scores of rate parameters or expression changes

# Examples
```julia
# Classify Z-scores
classify_Zscore(2.5, 2.0)   # Returns: "Up"
classify_Zscore(-1.5, 2.0)  # Returns: "Null"
classify_Zscore(-3.0, 2.0)  # Returns: "Down"
```
"""
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


"""
    make_dataframe_transient(folder::String, winners::String="")

Create a DataFrame from transient rate files in a folder.

# Arguments
- `folder::String`: Path to folder containing transient rate files
- `winners::String=""`: Optional path to winners file

# Returns
- `DataFrame`: DataFrame containing transient rate data

# Notes
- Searches for files containing "rate" in the filename
- Parses time information from filename (format: T{time}_...)
- Extracts rate parameters: on, off, eject, decay, yield
- Combines data from multiple time points
- Adds condition and time columns
- Optionally adds winner information if winners file is provided
- Converts time from minutes to hours (divides by 60)

# Examples
```julia
# Create transient DataFrame without winners
df = make_dataframe_transient("transient_results/")

# Create transient DataFrame with winners
df = make_dataframe_transient("transient_results/", "winners.txt")

# The resulting DataFrame contains:
# - Gene, on, off, eject, decay, yield, cond, time
# - winner (if winners file provided)
```
"""

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

"""
    filter_gene(measurefile, measure, threshold)

Filter genes based on a measure value and threshold.

# Arguments
- `measurefile::String`: Path to measure file
- `measure::String`: Name of the measure column
- `threshold::Float64`: Threshold value for filtering

# Returns
- `Vector{String}`: Vector of gene names that meet the criteria

# Notes
- Reads measure file and finds the specified measure column
- Returns genes where measure value > threshold OR is NaN
- Assumes first column contains gene names
- Prints the number of genes found for debugging
- Useful for identifying problematic genes in model fitting

# Examples
```julia
# Filter genes with high deviance
high_deviance_genes = filter_gene("measures.csv", "Deviance", 100.0)

# Filter genes with high Rhat values
high_rhat_genes = filter_gene("measures.csv", "Rhat", 1.1)
```
"""
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

"""
    filter_gene_nan(measurefile, measure)

Filter genes that have NaN values for a specific measure.

# Arguments
- `measurefile::String`: Path to measure file
- `measure::String`: Name of the measure column

# Returns
- `Vector{String}`: Vector of gene names with NaN values

# Notes
- Reads measure file and finds the specified measure column
- Returns genes where the measure value is NaN
- Assumes first column contains gene names
- Prints the number of genes found for debugging
- Useful for identifying genes with missing or failed model fits

# Examples
```julia
# Find genes with NaN deviance values
nan_deviance_genes = filter_gene_nan("measures.csv", "Deviance")

# Find genes with NaN AIC values
nan_aic_genes = filter_gene_nan("measures.csv", "AIC")
```
"""
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

"""
    ratestats(gene, G, folder, cond)

Get rate parameters and covariance matrix for a specific gene.

# Arguments
- `gene::String`: Gene name
- `G`: Model identifier
- `folder::String`: Path to folder containing result files
- `cond::String`: Condition identifier

# Returns
- `Tuple{Vector, Matrix}`: Tuple containing (rates, covariance_matrix)

# Notes
- Reads rate parameters from rate file
- Reads covariance matrix from parameter statistics file
- Assumes files follow standard naming conventions
- Uses getfile and getratefile functions to locate files
- Returns both rate parameters and their uncertainties

# Examples
```julia
# Get rate statistics for a gene
rates, cov = ratestats("MYC", "2", "results/", "control")

# rates contains the fitted parameters
# cov contains the covariance matrix for uncertainty analysis
```
"""
function ratestats(gene, G, folder, cond)
    filestats = joinpath(folder, getfile("param-stats", gene, G, folder, cond)[1])
    filerates = joinpath(folder, getratefile(gene, G, folder, cond)[1])
    rates = readrates(filerates)
    cov = read_covparam(filestats)
    # mu = readmean(filestats)
    # sd = readsd(filestats)
    return rates, cov
end

"""
    meanofftime(r::Vector, n::Int, method)

Calculate mean off time for a gene expression model.

# Arguments
- `r::Vector`: Rate parameters vector
- `n::Int`: Number of states
- `method`: Method for calculating off time

# Returns
- `Float64`: Mean off time

# Notes
- For single-state models (n=1): returns 1/r[1]
- For multi-state models: calculates sum of (1 - off time probabilities)
- Uses offtime function to get off time probabilities
- Mean off time represents average time spent in off states
- Important metric for understanding gene expression dynamics

# Examples
```julia
# Calculate mean off time for a two-state model
r = [0.1, 0.2, 0.5, 0.3]  # rate parameters
mean_off = meanofftime(r, 2, "method")

# Returns the average time the gene spends in off states
```
"""
function meanofftime(r::Vector, n::Int, method)
    if n == 1
        return 1 / r[1]
    else
        return sum(1 .- offtime(r, n, method))
    end
end

"""
    offtime(r::Vector, n::Int, method)

Calculate off time probabilities for a gene expression model.

# Arguments
- `r::Vector`: Rate parameters vector
- `n::Int`: Number of states
- `method`: Method for calculating off time

# Returns
- `Vector{Float64}`: Off time probabilities

# Notes
- Constructs gene state transition matrix using mat_G_DT
- Performs eigendecomposition to find eigenvalues
- Uses minimum eigenvalue to determine time range for CDF calculation
- Calculates off time cumulative distribution function
- Returns probabilities of being in off states at different times
- Important for understanding gene switching dynamics

# Examples
```julia
# Calculate off time probabilities for a two-state model
r = [0.1, 0.2, 0.5, 0.3]  # rate parameters
off_probs = offtime(r, 2, "method")

# Returns vector of off time probabilities
```
"""
function offtime(r::Vector, n::Int, method)
    _, _, TI = mat_G_DT(r, n)
    vals, _ = eig_decompose(TI)
    minval = min(minimum(abs.(vals[vals.!=0])), 0.2)
    offtimeCDF(collect(1.0:5/minval), r, n, TI, method)
end

"""
    offtime(gene::String, infile, n, method, root)

Calculate off time probabilities for a specific gene from a rate file.

# Arguments
- `gene::String`: Gene name
- `infile::String`: Path to rate file
- `n::Int`: Number of states
- `method`: Method for calculating off time
- `root::String`: Root path for data

# Returns
- `Vector{Float64}`: Off time probabilities

# Notes
- Reads rate parameters for the specified gene from the rate file
- Extracts rate parameters (excluding first and last columns)
- Delegates to offtime(r::Vector, n::Int, method) for calculation
- Assumes rate file has gene names in first column
- Assumes rate parameters are in columns 2 to end-1

# Examples
```julia
# Calculate off time for a gene from rate file
off_probs = offtime("MYC", "rates.txt", 2, "method", "data/")

# Returns off time probabilities for the MYC gene
```
"""
function offtime(gene::String, infile, n, method, root)
    contents, head = readdlm(infile, ',', header=true)
    r = float.(contents[gene.==contents[:, 1], 2:end-1])[1, :]
    offtime(r, n, method)
end

"""
    join_files(file1::String, file2::String, outfile::String, addlabel::Bool=true)

Join two rate files into a single output file.

# Arguments
- `file1::String`: Path to first rate file (G=2 model)
- `file2::String`: Path to second rate file (G=3 model)
- `outfile::String`: Path to output file
- `addlabel::Bool=true`: Whether to add model labels to column names

# Returns
- `Nothing`: Writes joined data to outfile

# Notes
- Joins files based on matching gene names in first column
- Adds "_G2" and "_G3" suffixes to column names if addlabel=true
- Preserves gene names in first column
- Only includes rows where gene names match between files
- Assumes both files have the same gene names in same order
- Useful for comparing results from different model types

# Examples
```julia
# Join G=2 and G=3 model results
join_files("rates_G2.txt", "rates_G3.txt", "combined_rates.txt")

# Join without model labels
join_files("rates_G2.txt", "rates_G3.txt", "combined_rates.txt", false)
```
"""

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

# function sample_non1_genes(infile, n)
#     contents, head = readdlm(infile, ',', header=true)
#     list = Array{String,1}(undef, 0)
#     for c in eachrow(contents)
#         if c[5] != 1
#             push!(list, c[1])
#         end
#     end
#     a = StatsBase.sample(list, n, replace=false)
# end

"""
    add_best_burst(filein, fileout, filemodel2, filemodel3)

Add burst size information to a file based on best model selection.

# Arguments
- `filein::String`: Path to input file containing model results
- `fileout::String`: Path to output file
- `filemodel2::String`: Path to G=2 model burst file
- `filemodel3::String`: Path to G=3 model burst file

# Returns
- `Nothing`: Writes augmented data to outfile

# Notes
- Reads model results and determines best model for each gene
- Adds burst size information based on the best model (G=2 or G=3)
- Assumes last column contains model identifier
- Adds "mean off period" and "burst size" columns
- Matches genes between input file and burst files
- Useful for adding burst statistics to model comparison results

# Examples
```julia
# Add burst information to model comparison results
add_best_burst("model_comparison.txt", "with_bursts.txt", "burst_G2.txt", "burst_G3.txt")
```
"""
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

"""
    add_best_occupancy(filein, fileout, filemodel2, filemodel3)

Add occupancy information to a file based on best model selection.

# Arguments
- `filein::String`: Path to input file containing model results
- `fileout::String`: Path to output file
- `filemodel2::String`: Path to G=2 model occupancy file
- `filemodel3::String`: Path to G=3 model occupancy file

# Returns
- `Nothing`: Writes augmented data to outfile

# Notes
- Reads model results and determines best model for each gene
- Adds occupancy information based on the best model (G=2 or G=3)
- Assumes third-to-last column contains model identifier
- Adds "Off -2", "Off -1", "On State" columns
- For G=2 models: adds "NA" for "Off -2" since G=2 only has one off state
- For G=3 models: adds all three occupancy values
- Matches genes between input file and occupancy files

# Examples
```julia
# Add occupancy information to model comparison results
add_best_occupancy("model_comparison.txt", "with_occupancy.txt", "occupancy_G2.txt", "occupancy_G3.txt")
```
"""
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

# function prune_file(list, file, outfile, header=true)
#     contents, head = readdlm(file, ',', header=header)
#     f = open(outfile, "w")
#     for c in eachrow(contents)
#         if c[1] in list
#             writedlm(f, [c], ',')
#         end
#     end
#     close(f)
# end

"""
    replace_yield(G, folder1, folder2, cond1, cond2, outfolder)

Replace yield parameters from one condition with those from another condition.

# Arguments
- `G`: Model identifier (converted to string if numeric)
- `folder1::String`: Path to source folder
- `folder2::String`: Path to target folder
- `cond1::String`: Source condition identifier
- `cond2::String`: Target condition identifier
- `outfolder::String`: Path to output folder

# Returns
- `Nothing`: Creates output files in outfolder

# Notes
- Reads rate files from both conditions
- Replaces the last parameter (yield) in target files with yield from source files
- Creates output folder if it doesn't exist
- Matches files by gene names between conditions
- Preserves all other parameters from target condition
- Useful for transferring yield estimates between experimental conditions

# Examples
```julia
# Replace yield from control condition into treatment condition
replace_yield("2", "control_results/", "treatment_results/", "control", "treatment", "modified_results/")
```
"""
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

"""
    getfile(type, gene::String, G::String, folder, cond)

Find files of a specific type for a gene, model, and condition.

# Arguments
- `type::String`: File type prefix (e.g., "rate", "param-stats")
- `gene::String`: Gene name
- `G::String`: Model identifier
- `folder::String`: Path to folder to search
- `cond::String`: Condition identifier

# Returns
- `Vector{String}`: Vector of matching file names

# Notes
- Searches for files containing type prefix, gene name, model, and condition
- Filters files by multiple criteria in sequence
- Returns all files that match all criteria
- Assumes files follow standard naming conventions
- Useful for finding specific file types for a gene-model-condition combination

# Examples
```julia
# Find rate files for MYC gene, G=2 model, control condition
files = getfile("rate", "MYC", "2", "results/", "control")

# Find parameter statistics files
files = getfile("param-stats", "MYC", "2", "results/", "control")
```
"""
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

"""
    make_halflife(infile, outfile, col=4)

Create a halflife file from input data with two halflife columns.

# Arguments
- `infile::String`: Path to input file
- `outfile::String`: Path to output file
- `col::Int=4`: Column index for first halflife value

# Returns
- `Nothing`: Writes halflife data to outfile

# Notes
- Reads input file and extracts halflife values from columns col and col+1
- Removes asterisks (*) from gene names
- Calculates average halflife when both values are available
- Only includes genes where at least one halflife value > 0
- Creates output file with "Gene" and "Halflife" columns
- Useful for processing halflife data from experimental measurements

# Examples
```julia
# Create halflife file from input data
make_halflife("input_data.txt", "halflife_output.txt")

# Use different column for first halflife value
make_halflife("input_data.txt", "halflife_output.txt", 6)
```
"""
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

"""
    make_datafiles(infolder, outfolder, label)

Copy files from input folder to output folder, removing a label from filenames.

# Arguments
- `infolder::String`: Path to input folder
- `outfolder::String`: Path to output folder
- `label::String`: Label to remove from filenames

# Returns
- `Nothing`: Copies files to outfolder

# Notes
- Creates output folder if it doesn't exist
- Copies all files from input folder to output folder
- Removes the specified label from each filename
- Preserves file content, only modifies filename
- Useful for cleaning up filenames by removing prefixes or suffixes

# Examples
```julia
# Copy files and remove "temp_" prefix
make_datafiles("input/", "output/", "temp_")

# This would copy and rename:
# - "temp_data1.txt" → "data1.txt"
# - "temp_data2.txt" → "data2.txt"
```
"""
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

"""
    make_observation_dist(d, states, G, R, S, coupling=tuple)

Create observation distributions for trace analysis based on model type and coupling.

# Arguments
- `d`: Observation distribution data
- `states`: State information for the model
- `G`: Gene states (Int or Tuple)
- `R`: RNA states (Int or Tuple)
- `S`: Splice states (Int or Tuple)
- `coupling`: Coupling structure (default: tuple())

# Returns
- `states`: Processed state information
- `observations`: Vector of observation distributions

# Notes
- Handles both uncoupled (G isa Int) and coupled (G isa Tuple) models
- For uncoupled models: processes single observation distributions
- For coupled models: processes multiple observation distributions with unit state mapping
- Automatically detects distribution types and handles accordingly
"""
function make_observation_dist(d, states, G, R, S, coupling=tuple)
    observations = Vector[]
    if G isa Int
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

"""
    make_trace_datamodel(traces::Vector, interval::Float64, rin, transitions, G, R, S, insertstep, probfn=prob_Gaussian, noiseparams=4, splicetype="", state=true, hierarchical=false, coupling=tuple(), grid=nothing, zeromedian=false)

Create trace data model from raw traces and model parameters.

# Arguments
- `traces::Vector`: Vector of trace data
- `interval::Float64`: Time interval between measurements
- `rin`: Rate parameters
- `transitions`: Model transition structure
- `G, R, S, insertstep`: Model parameters
- `probfn`: Probability function (default: prob_Gaussian)
- `noiseparams`: Number of noise parameters (default: 4)
- `splicetype`: Splice type (default: "")
- `state`: Include state information (default: true)
- `hierarchical`: Hierarchical model flag (default: false)
- `coupling`: Coupling structure (default: tuple())
- `grid`: Grid parameters (default: nothing)
- `zeromedian`: Zero median flag (default: false)

# Returns
- `data`: TraceData object
- `model`: Loaded model object

# Notes
- Creates TraceData object from raw traces
- Handles both hierarchical and non-hierarchical models
- Supports coupled and uncoupled model types
"""
function make_trace_datamodel(traces::Vector, interval::Float64, rin, transitions, G, R, S, insertstep, probfn=prob_Gaussian, noiseparams=4, splicetype="", state=true, hierarchical=false, coupling=tuple(), grid=nothing, zeromedian=false)
    data = TraceData{String,String,Tuple}("", "", interval, (traces, [], 0.0, length(traces[1])))
    make_trace_datamodel(data, rin, transitions, G, R, S, insertstep, probfn, noiseparams, splicetype, state, hierarchical, coupling, grid, zeromedian)
end

"""
    make_trace_datamodel(data::TraceData, rin, transitions, G, R, S, insertstep, probfn=prob_Gaussian, noiseparams=4, splicetype="", state=true, hierarchical=false, coupling=tuple(), grid=nothing, zeromedian=false)

Create trace data model from existing TraceData object and model parameters.

# Arguments
- `data::TraceData`: Existing trace data object
- `rin`: Rate parameters
- `transitions`: Model transition structure
- `G, R, S, insertstep`: Model parameters
- `probfn`: Probability function (default: prob_Gaussian)
- `noiseparams`: Number of noise parameters (default: 4)
- `splicetype`: Splice type (default: "")
- `state`: Include state information (default: true)
- `hierarchical`: Hierarchical model flag (default: false)
- `coupling`: Coupling structure (default: tuple())
- `grid`: Grid parameters (default: nothing)
- `zeromedian`: Zero median flag (default: false)

# Returns
- `data`: Updated TraceData object
- `model`: Loaded model object

# Notes
- Uses existing TraceData object instead of creating new one
- Handles both hierarchical and non-hierarchical models
- Supports coupled and uncoupled model types
"""
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


"""
    make_traces_dataframe(ts, td, traces, G::Int, R::Int, S::Int, insertstep::Int, state::Bool, coupling)

Create DataFrame from trace data for uncoupled models (Int parameters).

# Arguments
- `ts`: State trajectories
- `td`: Trace data
- `traces`: Original trace data
- `G::Int`: Number of gene states
- `R::Int`: Number of RNA states
- `S::Int`: Number of splice states
- `insertstep::Int`: Insertion step
- `state::Bool`: Include state information
- `coupling`: Coupling structure

# Returns
- `DataFrame`: DataFrame containing trace data, model predictions, and optional state information

# Notes
- Handles uncoupled models with integer parameters
- Creates columns for data, model means, and standard deviations
- Optionally includes gene states, RNA states, and reporter information
- Handles missing values by padding with missing
"""
function make_traces_dataframe(ts, td, traces, G::Int, R::Int, S::Int, insertstep::Int, state::Bool, coupling)
    l = maximum(length.(traces))
    if !isempty(coupling)
        data = ["data_$i" => [traces[i][:, 2]; fill(missing, l - length(traces[i][:, 2]))] for i in eachindex(traces)]
    else
        data = ["data_$i" => [traces[i]; fill(missing, l - length(traces[i]))] for i in eachindex(traces)]
    end
    pred = ["model_mean_$i" => [mean.(td[i]); fill(missing, l - length(td[i]))] for i in eachindex(td)]
    predstd = ["model_std_$i" => [std.(td[i]); fill(missing, l - length(td[i]))] for i in eachindex(td)]
    cols = [data pred predstd]
    if state
        g, z, zdigits, r = inverse_state(ts, G, R, S, insertstep)
        gs = ["Gstate_$i" => [g[i]; fill(missing, l - length(g[i]))] for i in eachindex(g)]
        s = ["Rstate_$i" => [zdigits[i]; fill(missing, l - length(g[i]))] for i in eachindex(g)]
        zs = ["Reporters_$i" => [r[i]; fill(missing, l - length(g[i]))] for i in eachindex(g)]
        cols = hcat(cols, [gs s zs])
    end
    DataFrame(permutedims(cols, (2, 1))[:])
end

"""
    make_traces_dataframe(ts, tp, traces, G::Tuple, R::Tuple, S::Tuple, insertstep::Tuple, state::Bool, coupling)

Create DataFrame from trace data for coupled models (Tuple parameters).

# Arguments
- `ts`: State trajectories
- `tp`: Trace predictions
- `traces`: Original trace data
- `G::Tuple`: Gene states for each unit
- `R::Tuple`: RNA states for each unit
- `S::Tuple`: Splice states for each unit
- `insertstep::Tuple`: Insertion steps for each unit
- `state::Bool`: Include state information
- `coupling`: Coupling structure

# Returns
- `DataFrame`: DataFrame containing trace data, model predictions, and optional state information

# Notes
- Handles coupled models with tuple parameters
- Creates separate columns for each coupled unit
- Processes each unit in the coupling structure
- Optionally includes state information for each unit
- Handles missing values by padding with missing
"""
function make_traces_dataframe(ts, tp, traces, G::Tuple, R::Tuple, S::Tuple, insertstep::Tuple, state::Bool, coupling)
    l = maximum(size.(traces, 1))
    cols = Matrix(undef, length(traces), 0)
    for k in coupling[1]
        data = ["data$k" * "_$i" => [traces[i][:, k]; fill(missing, l - length(traces[i][:, k]))] for i in eachindex(traces)]
        pred = ["model_mean$k" * "_$i" => [[mean(t[k]) for t in tp[i]]; fill(missing, l - length(tp[i]))] for i in eachindex(tp)]
        predstd = ["std_mean$k" * "_$i" => [[std(t[k]) for t in tp[i]]; fill(missing, l - length(tp[i]))] for i in eachindex(tp)]
        cols = hcat(cols, [data pred predstd])
        if state
            index = [[s[k] for s in t] for t in ts]
            g, z, zdigits, r = inverse_state(index, G[k], R[k], S[k], insertstep[k])
            gs = ["Gstate$k" * "_$i" => [g[i]; fill(missing, l - length(g[i]))] for i in eachindex(g)]
            s = ["Rstate$k" * "_$i" => [zdigits[i]; fill(missing, l - length(g[i]))] for i in eachindex(g)]
            zs = ["Reporters$k" * "_$i" => [r[i]; fill(missing, l - length(g[i]))] for i in eachindex(g)]
            cols = hcat(cols, [gs s zs])
        end
    end
    DataFrame(permutedims(cols, (2, 1))[:])
end

"""
    make_traces_dataframe(traces, interval, rin, transitions, G, R, S, insertstep, probfn=prob_Gaussian, noiseparams=4, splicetype="", state=true, hierarchical=false, coupling=tuple(), grid=nothing, zeromedian=false)

Create DataFrame from raw traces and model parameters.

# Arguments
- `traces`: Vector of trace data
- `interval`: Time interval between measurements
- `rin`: Rate parameters
- `transitions`: Model transition structure
- `G, R, S, insertstep`: Model parameters
- `probfn`: Probability function (default: prob_Gaussian)
- `noiseparams`: Number of noise parameters (default: 4)
- `splicetype`: Splice type (default: "")
- `state`: Include state information (default: true)
- `hierarchical`: Hierarchical model flag (default: false)
- `coupling`: Coupling structure (default: tuple())
- `grid`: Grid parameters (default: nothing)
- `zeromedian`: Zero median flag (default: false)

# Returns
- `DataFrame`: DataFrame containing trace data and model predictions

# Notes
- Creates TraceData object from raw traces
- Delegates to the TraceData version of make_traces_dataframe
- Handles both coupled and uncoupled models automatically
"""
function make_traces_dataframe(traces, interval, rin, transitions, G, R, S, insertstep, probfn=prob_Gaussian, noiseparams=4, splicetype="", state=true, hierarchical=false, coupling=tuple(), grid=nothing, zeromedian=false)
    data = TraceData{String,String,Tuple}("", "", interval, (traces, [], 0.0, length(traces[1])))
    make_traces_dataframe(data, rin, transitions, G, R, S, insertstep, probfn, noiseparams, splicetype, state, hierarchical, coupling, grid, zeromedian)
end

"""
    make_traces_dataframe(data::AbstractTraceData, rin, transitions, G, R, S, insertstep, probfn=prob_Gaussian, noiseparams=4, splicetype="", state=true, hierarchical=false, coupling=tuple(), grid=nothing, zeromedian=false)

Create DataFrame from TraceData object and model parameters.

# Arguments
- `data::AbstractTraceData`: Trace data object
- `rin`: Rate parameters
- `transitions`: Model transition structure
- `G, R, S, insertstep`: Model parameters
- `probfn`: Probability function (default: prob_Gaussian)
- `noiseparams`: Number of noise parameters (default: 4)
- `splicetype`: Splice type (default: "")
- `state`: Include state information (default: true)
- `hierarchical`: Hierarchical model flag (default: false)
- `coupling`: Coupling structure (default: tuple())
- `grid`: Grid parameters (default: nothing)
- `zeromedian`: Zero median flag (default: false)

# Returns
- `DataFrame`: DataFrame containing trace data and model predictions

# Notes
- Loads model and generates predictions
- Handles both hierarchical and non-hierarchical models
- Distinguishes between coupled (G isa Tuple) and uncoupled models
- Creates observation distributions and state information
- Delegates to appropriate make_traces_dataframe method based on parameter types
"""
function make_traces_dataframe(data::AbstractTraceData, rin, transitions, G, R, S, insertstep, probfn=prob_Gaussian, noiseparams=4, splicetype="", state=true, hierarchical=false, coupling=tuple(), grid=nothing, zeromedian=false)
    # data = TraceData{String,String,Tuple}("", "", interval, (traces, [], 0.0, length(traces[1])))
    if hierarchical
        h = (2, [8], ())
        method = (Tsit5(), true)
    else
        h = ()
        method = Tsit5()
    end
    if !isempty(coupling) && G isa Tuple
        model = load_model(data, rin, rin, [1, 2, 3], (), transitions, G, R, S, insertstep, splicetype, 1, 10.0, Int[], 1.0, 0.1, probfn, [ones(Int, noiseparams), ones(Int, noiseparams)], method, h, coupling, grid, zeromedian)
    else
        model = load_model(data, rin, rin, [1, 2, 3], (), transitions, G, R, S, insertstep, splicetype, 1, 10.0, Int[], 1.0, 0.1, probfn, ones(Int, noiseparams), method, h, coupling, grid, zeromedian)
    end
    ts, d = predict_trace(get_param(model), data, model)
    states, observations = make_observation_dist(d, ts, G, R, S, coupling)
    make_traces_dataframe(states, observations, data.trace[1], G, R, S, insertstep, state, coupling)
end

function write_trace_dataframe(outfile, datapath, datacond, interval::Float64, r::Vector, transitions, G, R, S, insertstep, start=1, stop=-1, probfn=prob_Gaussian, noiseparams=4, splicetype=""; state=true, hierarchical=false, coupling=tuple(), grid=nothing, zeromedian=false, datacol=3)
    traceinfo = (interval, start, stop, 1.0, 0.0)
    data = load_data_trace(datapath, "", "", datacond, traceinfo, :trace, datacol, zeromedian)
    df = make_traces_dataframe(data, r, transitions, G, R, S, insertstep, probfn, noiseparams, splicetype, state, hierarchical, coupling, grid, zeromedian)
    CSV.write(outfile, df)
end


function write_trace_dataframe(file::String, datapath::String, interval::Float64, ratetype::String="median", start=1, stop=-1, probfn=prob_Gaussian, noiseparams=4, splicetype=""; hlabel="-h", state=true, grid=nothing, zeromedian=false, datacol=3)
    println(file)
    datacond, transitions, G, R, S, insertstep, hierarchical, coupling_field = parse_filename(file, hlabel=hlabel)
    coupling = make_coupling(coupling_field, G, R)
    if G isa Int && !isempty(coupling)
        datapath = joinpath(datapath, string(coupling_field[1]))
    end
    println(datapath)
    r = readrates(file, get_row(ratetype))
    out = replace(file, "rates" => "predictedtraces", ".txt" => ".csv")
    write_trace_dataframe(out, datapath, datacond, interval, r, transitions, G, R, S, insertstep, start, stop, probfn, noiseparams, splicetype, state=state, hierarchical=hierarchical, coupling=coupling, grid=grid, zeromedian=zeromedian, datacol=datacol)
end

function write_traces(folder::String, datapath::String, interval::Float64, ratetype::String="median", start=1, stop=-1, probfn=prob_Gaussian, noiseparams=4, splicetype=""; hlabel="-h", state=true, grid=nothing, zeromedian=false, datacol=3)
    # Collect all files to process first
    files_to_process = String[]
    for (root, dirs, files) in walkdir(folder)
        for f in files
            if occursin("rates", f)
                push!(files_to_process, joinpath(root, f))
            end
        end
    end
    
    # Process files in parallel
    Threads.@threads for file_path in files_to_process
        write_trace_dataframe(file_path, datapath, interval, ratetype, start, stop, probfn, noiseparams, splicetype, hlabel=hlabel, state=state, grid=grid, zeromedian=zeromedian, datacol=datacol)
    end
end



# """
#     write_trace_dataframe(file, datapath::String, datacond, interval, ratetype::String="median", start=1, stop=-1, probfn=prob_Gaussian, noiseparams=4, splicetype=""; hlabel="-h", state=true, coupling=tuple())

# TBW
# """
# function write_trace_dataframe(outfile, datapath, datacond, interval::Float64, r::Vector, transitions, G, R, S, insertstep, start=1, stop=-1, probfn=prob_Gaussian, noiseparams=4, splicetype=""; state=true, hierarchical=false, coupling=tuple(), grid=nothing, zeromedian=false)
#     traces = read_tracefiles(datapath, datacond, start, stop)
#     traces, maxmedians = zero_median(traces, zeromedian)
#     df = make_traces_dataframe(traces, interval, r, transitions, G, R, S, insertstep, probfn, noiseparams, splicetype, state, hierarchical, coupling, grid, zeromedian)
#     CSV.write(outfile, df)
# end
# function write_trace_dataframe(file, datapath::String, datacond, interval, ratetype::String="median", start=1, stop=-1, probfn=prob_Gaussian, noiseparams=4, splicetype=""; hlabel="-h", state=true, coupling=tuple(), grid=nothing, zeromedian=false, datacol=3)
#     println(file)
#     filename = basename(file)
#     occursin(hlabel, filename) ? hierarchical = true : hierarchical = false
#     parts = fields(filename)
#     G, R, S, insertstep = decompose_model(parts.model)
#     transitions = get_transitions(G, parts.label)
#     r = readrates(file, get_row(ratetype))
#     out = replace(file, "rates" => "predictedtraces", ".txt" => ".csv")
#     write_trace_dataframe(out, datapath, datacond, interval, r, transitions, G, R, S, insertstep, start, stop, probfn, noiseparams, splicetype, state=state, hierarchical=hierarchical, coupling=coupling, grid=grid, zeromedian=zeromedian, datacol=datacol)
# end
# """
#     write_traces(folder, datapath::String, datacond, interval, ratetype::String="median", start=1, stop=-1, probfn=prob_Gaussian, noiseparams=4, splicetype=""; hlabel="-h", state=true, coupling=tuple())

# """
# function write_traces(folder, datapath::String, datacond, interval, ratetype::String="median", start=1, stop=-1, probfn=prob_Gaussian, noiseparams=4, splicetype=""; hlabel="-h", state=true, coupling=tuple(), grid=nothing, zeromedian=false, datacol=3)
#     for (root, dirs, files) in walkdir(folder)
#         for f in files
#             # if occursin("rates", f) && occursin(datacond, f) #&& ((!exclude_label && occursin(hlabel, f)) || exclude_label && !occursin(hlabel, f))
#             if occursin("rates", f) && (occursin("tracejoint", f) || (isempty(coupling) && occursin(datacond, f)))
#                 write_trace_dataframe(joinpath(root, f), datapath, datacond, interval, ratetype, start, stop, probfn, noiseparams, splicetype, hlabel=hlabel, state=state, coupling=coupling, grid=grid, zeromedian=zeromedian)
#             end
#         end
#     end
# end
# """
#     write_traces_folder(folder, datafolder::Vector, datacond, interval, ratetype::String="median", start=1, stop=-1, probfn=prob_Gaussian, noiseparams=4, splicetype=""; state=true, coupling=tuple())

# TBW
# """
# function write_traces_coupling(folder, datapath, datacond, interval, G=(3, 3), R=(3, 3), sources=1:3, targets=1:5, ratetype::String="median", start=1, stop=-1, probfn=prob_Gaussian, noiseparams=4, splicetype=""; hlabel="-h", state=true, pattern="gene", zeromedian=true, datacol=3)
#     for (root, dirs, files) in walkdir(folder)
#         for f in files
#             for target in targets
#                 for source in sources
#                     # if occursin("rates", f) && occursin(datacond, f) #&& ((!exclude_label && occursin(hlabel, f)) || exclude_label && !occursin(hlabel, f))
#                     if occursin("rates", f) && occursin("$pattern$source$target", f)
#                         write_trace_dataframe(joinpath(root, f), datapath, datacond, interval, ratetype, start, stop, probfn, noiseparams, splicetype, hlabel=hlabel, state=state, coupling=((1, 2), (tuple(), tuple(1)), (source, 0), (0, target), 1), zeromedian=zeromedian, datacol=datacol)
#                     end
#                 end
#                 if occursin("rates", f) && occursin("R$target", f)
#                     coupling = ((1, 2), (tuple(), tuple(1)), (collect(G[1]+1:G[1]+R[1]), 0), (0, target), 1)
#                     write_trace_dataframe(joinpath(root, f), datapath, datacond, interval, ratetype, start, stop, probfn, noiseparams, splicetype, hlabel=hlabel, state=state, coupling=coupling, zeromedian=zeromedian, datacol=datacol)
#                 end
#             end
#         end
#     end
# end



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
# """
#     write_traces_coupling_spawn(folder, datapath, datacond, interval, G=(3, 3), R=(3, 3),
#                                sources=1:3, targets=1:5, ratetype::String="median",
#                                start=1, stop=-1, probfn=prob_Gaussian, noiseparams=4,
#                                splicetype=""; hlabel="-h", state=true, pattern="gene")

# Parallelized version of write_traces_coupling using Julia's task-based parallelism with @spawn.
# This provides more dynamic scheduling which can be beneficial for workloads with varying execution times.
# """
# function write_traces_coupling_spawn(folder, datapath, datacond, interval, G=(3, 3), R=(3, 3),
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
#                             file_path = joinpath(root, f)
#                             coupling = ((1, 2), (tuple(), tuple(1)), (source, 0), (0, target), 1)

#                             # Create a task for this file
#                             t = @spawn begin
#                                 try
#                                     write_trace_dataframe(file_path, datapath, datacond, interval, ratetype,
#                                         start, stop, probfn, noiseparams, splicetype,
#                                         hlabel=hlabel, state=state, coupling=coupling)
#                                     @info "Successfully processed: $(basename(file_path))"
#                                 catch e
#                                     @error "Error processing file: $(basename(file_path))" exception = (e, catch_backtrace())
#                                 end
#                             end

#                             push!(tasks, t)
#                         end
#                     end

#                     # Tasks for R$target files
#                     if occursin("R$target", f)
#                         file_path = joinpath(root, f)
#                         coupling = ((1, 2), (tuple(), tuple(1)), (collect(G[1]+1:G[1]+R[1]), 0), (0, target), 1)

#                         # Create a task for this file
#                         t = @spawn begin
#                             try
#                                 write_trace_dataframe(file_path, datapath, datacond, interval, ratetype,
#                                     start, stop, probfn, noiseparams, splicetype,
#                                     hlabel=hlabel, state=state, coupling=coupling)
#                                 @info "Successfully processed: $(basename(file_path))"
#                             catch e
#                                 @error "Error processing file: $(basename(file_path))" exception = (e, catch_backtrace())
#                             end
#                         end

#                         push!(tasks, t)
#                     end
#                 end
#             end
#         end
#     end

#     @info "Scheduled $(length(tasks)) files for processing using $(Threads.nthreads()) threads"

#     # Wait for all tasks to complete
#     for t in tasks
#         fetch(t)
#     end

#     @info "All files processed"
# end

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



"""
    extract_source_target(pattern::String, filepath::String)

Extract source and target information from a filepath based on a pattern.

# Arguments
- `pattern::String`: Pattern to search for in filename
- `filepath::String`: Path to file

# Returns
- `Union{Tuple{Union{String,Int},Int}, Nothing}`: (source, target) or nothing if not found

# Notes
- Searches for pattern followed by digits in filename
- Handles both numeric sources and "R" (RNA) sources
- Returns source as integer or "R" string, target as integer
- Uses regex matching to find pattern## format
- Returns nothing if pattern is not found
- Useful for parsing coupling information from filenames

# Examples
```julia
# Extract source and target from filename
source, target = extract_source_target("gene", "rates_gene12_condition.txt")
# Returns: (1, 2)

# Extract RNA target
source, target = extract_source_target("gene", "rates_geneR3_condition.txt")
# Returns: ("R", 3)
```
"""
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



"""
    write_residency_G_allgenes(fileout::String, filein::String, G, header)

Compute gene state residence probabilities for all genes in a rate file and write results to CSV.

# Arguments
- `fileout::String`: Path to the output file
- `filein::String`: Path to the input file
- `G::Int`: Number of gene states
- `header::String`: Header for the output file

# Returns
- `Nothing`

# Notes
- Reads all rate parameters from the input file
- Computes residence probabilities for each gene
- Writes results to the output file
"""
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

"""
    write_residency_G(file, ratetype="median", transitions=..., nnoiseparams=4)

Compute gene state residence probabilities from a rate file and write results to CSV.

This function calculates the steady-state probability of being in each gene state (residence
probabilities) based on fitted rate parameters. Residence probabilities represent the fraction
of time the gene spends in each state at equilibrium, which is a key metric for understanding
gene expression dynamics.

# Arguments
- `file::String`: Path to the rate file (typically a `.txt` file containing fitted rate parameters).
  The filename should follow the standard naming convention so that model parameters (G, R, S, insertstep)
  can be extracted via `decompose_model`.
- `ratetype::String`: Which row of rates to use from the rate file (default: `"median"`). Other options
  include `"ml"` (maximum likelihood), `"mean"`, etc. This selects which parameter estimate to use when
  multiple estimates are stored in the file.
- `transitions::Tuple`: Tuple of transition definitions for each unit. Only used for coupled models
  (when `G` is a tuple). Default: `(([1, 2], [2, 1], [2, 3], [3, 2]), ([1, 2], [2, 1], [2, 3], [3, 2]))`
  for a 3-state model. For 2-state models (G=(2,2)), this is automatically set to `(([1, 2], [2, 1]), ([1, 2], [2, 1]))`.
- `nnoiseparams::Int`: Number of noise parameters in the rate vector (default: 4). Used to correctly
  parse rate parameters for coupled models.

# Output File

Creates a CSV file with the same name as the input file but with:
- `"rates"` replaced by `"residencyG"`
- `.txt` extension replaced by `.csv`

For example: `rates_trace-TEST-nstate_testtrace_test_3221_1.txt` → `residencyG_trace-TEST-nstate_testtrace_test_3221_1.csv`

## Output Format

The CSV file contains a single row with columns `G1`, `G2`, ..., `GN` where `N` is the number of gene states.
Each value represents the steady-state probability of being in that gene state.

**Example output:**
```csv
G1,G2,G3
0.300514,0.332833,0.366653
```

# Algorithm

For a single-gene model (G is an Int):
1. Reads rate parameters from the file using `readrates(file, get_row(ratetype))`.
2. Extracts model structure (G, R, S, insertstep) from the filename using `decompose_model`.
3. Computes residence probabilities using `residenceprob_G`, which solves the steady-state
   equations: `P(G_i) = P(G_{i-1}) × (rate_{i-1→i} / rate_{i→i-1})` for i > 1, with normalization.
4. Writes results to CSV with state labels as column headers.

For coupled models (G is a Tuple):
1. Determines the number of rate parameters per unit using `num_rates(transitions, R, S, insertstep) + nnoiseparams`.
2. Computes residence probabilities separately for each unit.
3. Writes results with columns labeled `G1_1`, `G2_1`, ... for unit 1 and `G1_2`, `G2_2`, ... for unit 2.

# Usage Example

```julia
# Compute residence probabilities for a single rate file
write_residency_G("results/trace-test/rates_trace-TEST-nstate_testtrace_test_3221_1.txt")

# Read and display results
using CSV, DataFrames
df = CSV.read("results/trace-test/residencyG_trace-TEST-nstate_testtrace_test_3221_1.csv", DataFrame)
println(df)
# 1×3 DataFrame
#  Row │ G1        G2        G3
#      │ Float64   Float64   Float64
# ─────┼──────────────────────────────
#    1 │ 0.300514  0.332833  0.366653
```

# Notes

- **Steady-State Assumption**: Residence probabilities assume the system is at equilibrium. For
  time-dependent analysis, use simulation-based methods.

- **Rate Ordering**: The function assumes rates are ordered as: G transitions, R transitions, S transitions,
  decay rates, and noise parameters. This ordering must match the model structure.

- **Coupled Models**: For coupled models (two units), the function automatically handles the rate
  parameter parsing based on the number of states and transitions for each unit.

- **Model Extraction**: The function relies on the filename containing model information in a parseable
  format. If `decompose_model` fails, the function may error.

# See Also
- `write_residency_G_folder`: Batch process all rate files in a folder
- `residenceprob_G`: Core function that computes residence probabilities from rates
- `residenceprob_G_dataframe`: Formats residence probabilities as a DataFrame
- `decompose_model`: Extracts model parameters from filename
"""
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

Create histogram from trace data for a specific condition.

# Arguments
- `datapath::String`: Path to data directory
- `datacond::String`: Data condition identifier
- `start::Int=1`: Starting index for trace data
- `stop::Int=-1`: Ending index for trace data (-1 for all)

# Returns
- `Tuple{Vector, Any}`: (trace_data, histogram)

# Notes
- Reads trace files for the specified condition
- Combines all trace data into a single vector
- Creates normalized histogram from combined data
- Uses read_tracefiles to load data
- Uses histogram function with normalize=true
- Useful for analyzing overall distribution of trace values

# Examples
```julia
# Create histogram from trace data
trace_data, hist = make_trace_histogram("data/", "control")

# hist contains normalized histogram of all trace values
```
"""
function make_trace_histogram(datapath, datacond, start=1, stop=-1)
    traces = read_tracefiles(datapath, datacond, start, stop)
    ft = reduce(vcat, traces)
    h = histogram(ft, normalize=true)
    return ft, h
end

"""
    reset_noise!(r, transitions, G::Tuple, R, S, insertstep, num_noiseparams)

Reset noise parameters in a rate vector to default values.

# Arguments
- `r::Vector`: Rate parameter vector (modified in place)
- `transitions`: Model transition structure
- `G::Tuple`: Gene states for each unit
- `R`: RNA states for each unit
- `S`: Splice states for each unit
- `insertstep`: Insertion steps for each unit
- `num_noiseparams::Int`: Number of noise parameters per unit

# Returns
- `Nothing`: Modifies rate vector in place

# Notes
- Resets noise parameters to default values: [20, 0.1, 80, 0.1]
- Processes each unit in the coupled model
- Calculates position of noise parameters based on model structure
- Assumes noise parameters follow rate parameters for each unit
- Useful for reinitializing noise parameters during model fitting

# Examples
```julia
# Reset noise parameters in coupled model
reset_noise!(r, transitions, (3, 3), (3, 3), (1, 0), (1, 1), 4)

# Noise parameters are reset to default values
```
"""
function reset_noise!(r, transitions, G::Tuple, R, S, insertstep, num_noiseparams)
    n = 0
    for i in eachindex(G)
        n += num_rates(transitions[i], R[i], S[i], insertstep[i])
        r[n+1:n+num_noiseparams] .= [20, 0.1, 80, 0.1]
        n += num_noiseparams
    end
end

"""
    delete_noise!(r, transitions, G::Tuple, R, S, insertstep, num_noiseparams)

Remove noise parameters from a rate vector.

# Arguments
- `r::Vector`: Rate parameter vector (modified in place)
- `transitions`: Model transition structure
- `G::Tuple`: Gene states for each unit
- `R`: RNA states for each unit
- `S`: Splice states for each unit
- `insertstep`: Insertion steps for each unit
- `num_noiseparams::Int`: Number of noise parameters per unit

# Returns
- `Nothing`: Modifies rate vector in place

# Notes
- Removes noise parameters from the rate vector
- Processes each unit in the coupled model
- Calculates position of noise parameters based on model structure
- Assumes noise parameters follow rate parameters for each unit
- Useful for creating rate vectors without noise parameters

# Examples
```julia
# Remove noise parameters from coupled model
delete_noise!(r, transitions, (3, 3), (3, 3), (1, 0), (1, 1), 4)

# Rate vector now contains only rate parameters
```
"""
function delete_noise!(r, transitions, G::Tuple, R, S, insertstep, num_noiseparams)
    n = 0
    for i in eachindex(G)
        n += num_rates(transitions[i], R[i], S[i], insertstep[i])
        deleteat!(r, n+1:n+num_noiseparams)
        # n += num_noiseparams
    end
end


"""
    sample_non1_genes(infile, n)

Sample genes that have non-unity values in column 5.

# Arguments
- `infile::String`: Path to input file
- `n::Int`: Number of genes to sample

# Returns
- `Vector{String}`: Vector of sampled gene names

# Notes
- Reads the input file and filters genes where column 5 != 1
- Assumes gene names are in the first column
- Uses StatsBase.sample for random sampling without replacement
- Useful for selecting a subset of genes for analysis
- Returns empty vector if no genes meet criteria

# Examples
```julia
# Sample 10 genes with non-unity values in column 5
genes = sample_non1_genes("data.txt", 10)

# Returns a random sample of gene names
```
"""
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

"""
    prune_file(list, file, outfile, header=true)

Create a new file containing only genes from a specified list.

# Arguments
- `list::Vector{String}`: List of gene names to keep
- `file::String`: Path to input file
- `outfile::String`: Path to output file
- `header::Bool=true`: Whether the input file has a header

# Returns
- `Nothing`: Writes filtered data to outfile

# Notes
- Reads the input file and filters rows based on gene names
- Assumes gene names are in the first column
- Only keeps rows where the gene name is in the provided list
- Preserves header if header=true
- Useful for creating subset files for specific gene sets
- Writes results to the specified output file

# Examples
```julia
# Keep only specific genes
genes_to_keep = ["MYC", "FOS", "JUN"]
prune_file(genes_to_keep, "all_genes.txt", "subset_genes.txt")

# Process file without header
prune_file(genes_to_keep, "data.txt", "subset.txt", header=false)
```
"""
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

"""
    split_conditions(cond::AbstractString, multicond::Bool)

Split condition string based on multicond flag.

# Arguments
- `cond::AbstractString`: Condition string to split
- `multicond::Bool`: If true, split on "-" delimiter; if false, return as single-element vector

# Returns
- `Vector{String}`: Vector of condition strings

# Examples
```julia
split_conditions("ctrl-treated", true)   # Returns: ["ctrl", "treated"]
split_conditions("ctrl", false)          # Returns: ["ctrl"]
```
"""
split_conditions(cond::AbstractString, multicond::Bool) = multicond ? split(cond, "-") : [cond]


################################################################################
# simulate_trials function
################################################################################


"""
    simulate_trials(r, transitions, G, R, S, insertstep, coupling, ntrials; trial_time=720.0, lags=collect(0:60), probfn=prob_Gaussian, correlation_algorithm=StandardCorrelation(), warmupsteps=1000)

Simulate multiple trials and compute empirical correlation functions from simulated traces, comparing to theoretical predictions.

The empirical correlation functions are computed using the specified `correlation_algorithm` (CorrelationTrait), which determines
centering, normalization, multi-tau binning, and biased/unbiased divisors. The theory is raw uncentered unnormalized correlation
functions (E[xy]) from the HMM, which are then transformed to match the empirical CorrelationTrait for comparison.

This function integrates simulation with empirical covariance computation:
1. Computes theoretical cross-covariances from model parameters
2. Simulates multiple traces from the model
3. Computes empirical cross-correlation functions from simulated traces using `compute_correlation_functions`
4. Compares theoretical and empirical results

# Arguments
- `r::Vector{Float64}`: Rate parameters for the model
- `transitions::Tuple`: State transition structure
- `G, R, S, insertstep`: Model structure parameters
- `coupling::Tuple`: Coupling structure for coupled models
- `ntrials::Int`: Number of simulation trials to run
- `trial_time::Float64=720.0`: Length of each simulated trace (minutes)
- `lag::Int=60`: Maximum lag for cross-covariance computation
- `stride::Int=1`: Step size for lags (default: 1 for every lag)
- `probfn::Function=prob_Gaussian`: Probability function for noise

# Returns
NamedTuple with:
- `cc_theory`, `cc_mean`: Theoretical and empirical intensity cross-covariance
- `ccON_theory`, `ccON_mean`: Theoretical and empirical ON state cross-covariance
- `ac1_mean`, `ac2_mean`: Empirical autocovariances
- `ac1ON_mean`, `ac2ON_mean`: Empirical ON state autocovariances
- `linf_norm`, `l2_norm`: Norms comparing theory vs empirical
- `mean1`, `mean2`, `m1`, `m2`: Mean values
- `v1`, `v2`, `v1_empirical`, `v2_empirical`: Variance values
- `lags`: Time lags used
- Bootstrap confidence intervals if available (`cc_se`, `cc_lower`, etc.)

# Notes
- Uses `compute_correlation_functions` internally to compute empirical cross-correlations from simulated traces
- Performs bootstrap resampling for error estimation
- Compares theoretical predictions (from `covariance_functions`) with empirical results from simulations
- Useful for validating model predictions and assessing finite-sample effects

# Example
```julia
r = [0.1, 0.2, 0.5, 1.0]  # Rate parameters
transitions = (([1, 2], [2, 1]), ([1, 2], [2, 1]))
G, R, S = (3, 3), (3, 3), (0, 0)
coupling = tuple()  # No coupling

result = simulate_trials(r, transitions, G, R, S, (1, 1), coupling, 100, 
                         trial_time=720.0, lag=60)
println("L2 norm: ", result.l2_norm)
```
"""
function simulate_trials(r::Vector, transitions::Tuple, G, R, S, insertstep, coupling, ntrials, trial_time=720.0, lags=collect(0:60); probfn=prob_Gaussian, offset::Float64=0.0, correlation_algorithm=StandardCorrelation(), warmupsteps::Int=1000)
    # If lags is a single integer, treat it as max lag (backward compatibility)
    if lags isa Integer
        positive_lags = collect(0:lags)
    else
        # Extract positive lags for covariance_functions (it expects positive lags)
        lags_vec = collect(lags)
        positive_lags = unique(sort([abs(τ) for τ in lags_vec if τ >= 0]))
    end
    # Get both intensity and ON/OFF cross-covariances and autocovariances
    # Use same offset as in covariance_functions to ensure theory and empirical match
    # Note: interval is now inferred from lags, so we don't pass it
    full_lags, cc, ac1, ac2, m1, m2, v1, v2, ccON, ac1ON, ac2ON, m1ON, m2ON, v1ON, v2ON, ccReporters, ac1Reporters, ac2Reporters, m1Reporters, m2Reporters, v1Reporters, v2Reporters = correlation_functions(r, transitions, G, R, S, insertstep, probfn, coupling, positive_lags, offset=offset)
    # Pack into NamedTuple
    theory = (
        ac1=ac1, ac2=ac2, cc=cc,
        m1=m1, m2=m2, v1=v1, v2=v2,
        ccON=ccON, m1ON=m1ON, m2ON=m2ON, ac1ON=ac1ON, ac2ON=ac2ON,
        ccReporters=ccReporters, m1Reporters=m1Reporters, m2Reporters=m2Reporters,
        ac1Reporters=ac1Reporters, ac2Reporters=ac2Reporters,
        full_lags=full_lags
    )
    # Use ccON (ON/OFF cross-covariance) and ac1ON/ac2ON (ON/OFF autocovariances) for comparison
    # Use the provided lags (or convert if it was an integer)
    actual_lags = lags isa Integer ? full_lags : lags_vec
    simulate_trials(theory, r, transitions, G, R, S, insertstep, coupling, actual_lags, ntrials, trial_time, offset, correlation_algorithm, warmupsteps)
end

"""
    simulate_trials(theory, r, transitions, G, R, S, insertstep, coupling, lags_ac, lags; ntrials=1, trial_time=720.0, offset=0.0, correlation_algorithm=StandardCorrelation(), warmupsteps=1000)

Internal method that performs the actual simulation and comparison.

This is the lower-level method called by the main `simulate_trials` function.
It takes precomputed theoretical covariances and performs simulations.

# Arguments
- `theory::NamedTuple`: NamedTuple containing all theoretical correlation functions and means:
  - `ac1, ac2, cc`: Theoretical autocovariances and cross-covariance for intensity
  - `m1, m2, v1, v2`: Theoretical means and variances for intensity
  - `ccON, ac1ON, ac2ON`: Theoretical cross-covariance and autocovariances for ON states
  - `m1ON, m2ON`: Theoretical mean ON state probabilities
  - `ccReporters, ac1Reporters, ac2Reporters`: Theoretical cross-covariance and autocovariances for reporter counts
  - `m1Reporters, m2Reporters`: Theoretical mean reporter counts
  - `full_lags`: Full lag vector (symmetric)
- `r::Vector{Float64}`: Rate parameters (used for simulation)
- `transitions, G, R, S, insertstep, coupling`: Model structure parameters
- `lags_ac, lags::Vector{Int}`: Time lags for autocovariance and cross-covariance
- `ntrials::Int=1`: Number of simulation trials
- `trial_time::Float64=720.0`: Length of each simulated trace (minutes)
- `offset::Float64=0.0`: Offset for ON state computation
- `correlation_algorithm::CorrelationTrait`: Correlation algorithm to use for empirical computation
- `warmupsteps::Int=1000`: Number of warmup steps for simulation

# Returns
NamedTuple with theoretical and empirical results (see main `simulate_trials` docstring).
"""
function simulate_trials(theory, r::Vector, transitions::Tuple, G, R, S, insertstep, coupling, lags, ntrials=1, trial_time::Float64=720.0, offset::Float64=0.0, correlation_algorithm=StandardCorrelation(), warmupsteps::Int=1000)
    # Extract theory values from NamedTuple
    ac1 = theory.ac1
    ac2 = theory.ac2
    cc = theory.cc
    m1 = theory.m1
    m2 = theory.m2
    v1 = theory.v1
    v2 = theory.v2
    ccON = theory.ccON
    m1ON = theory.m1ON
    m2ON = theory.m2ON
    ac1ON = theory.ac1ON
    ac2ON = theory.ac2ON
    ccReporters = theory.ccReporters
    m1Reporters = theory.m1Reporters
    m2Reporters = theory.m2Reporters
    ac1Reporters = theory.ac1Reporters
    ac2Reporters = theory.ac2Reporters
    n_bootstrap = min(1000, max(50, ntrials * 10))

    intensity_traces = Vector{Matrix{Float64}}(undef, ntrials)
    reporter_traces = Vector{Matrix{Float64}}(undef, ntrials)
    on_traces = Vector{Matrix{Float64}}(undef, ntrials)

    Threads.@threads for i in 1:ntrials
        # Warmup is handled by simulator's warmupsteps parameter (converted from warmup_time)
        t = simulate_trace_vector(r, transitions, G, R, S, insertstep, coupling, 1.0, trial_time, 1, col=[2, 3], warmupsteps=warmupsteps)
        # Warmup is already handled by simulator, so use all the trace data
        trace_data = t[1]
        # trace_data columns: [intensity1, reporters1, intensity2, reporters2]
        intensity_traces[i] = hcat(Float64.(trace_data[:, 1]), Float64.(trace_data[:, 3]))
        reporter_traces[i] = hcat(Float64.(trace_data[:, 2]), Float64.(trace_data[:, 4]))
        # Apply same offset as theory: ON states = (binary > 0) + offset
        # This matches covariance_functions which computes: ON = float(num_per_state .> 0.0) .+ offset
        on_traces[i] = hcat(Float64.(trace_data[:, 2] .> 0.0) .+ offset, Float64.(trace_data[:, 4] .> 0.0) .+ offset)
    end

    # Compute empirical correlations using compute_covariance directly
    # lags should already be symmetric (full_lags from covariance_functions)
    # frame_interval = 1.0 because simulate_trace_vector uses interval=1.0 (one frame per minute)
    intensity_result = isempty(intensity_traces) ? nothing : compute_correlation_functions(intensity_traces, lags; correlation_algorithm=correlation_algorithm, bootstrap=true, n_bootstrap=n_bootstrap, frame_interval=1.0)
    on_result = compute_correlation_functions(on_traces, lags; correlation_algorithm=correlation_algorithm, bootstrap=true, n_bootstrap=n_bootstrap, frame_interval=1.0)
    reporter_result = compute_correlation_functions(reporter_traces, lags; correlation_algorithm=correlation_algorithm, bootstrap=true, n_bootstrap=n_bootstrap, frame_interval=1.0)

    # Extract results
    empirical = (
        cc=isnothing(intensity_result) ? Float64[] : intensity_result.cc,
        ac1=isnothing(intensity_result) ? Float64[] : intensity_result.ac1,
        ac2=isnothing(intensity_result) ? Float64[] : intensity_result.ac2,
        ccON=on_result.cc,
        ac1ON=on_result.ac1,
        ac2ON=on_result.ac2,
        ccReporters=reporter_result.cc,
        ac1Reporters=reporter_result.ac1,
        ac2Reporters=reporter_result.ac2,
        mON1=on_result.mean1,
        mON2=on_result.mean2,
        mR1=reporter_result.mean1,
        mR2=reporter_result.mean2,
        v1_empirical=isnothing(intensity_result) ? 0.0 : intensity_result.v1,
        v2_empirical=isnothing(intensity_result) ? 0.0 : intensity_result.v2,
        cc_lower=isnothing(intensity_result) ? nothing : intensity_result.cc_lower,
        cc_median=isnothing(intensity_result) ? nothing : intensity_result.cc_median,
        cc_upper=isnothing(intensity_result) ? nothing : intensity_result.cc_upper,
        cc_se=isnothing(intensity_result) ? nothing : intensity_result.cc_se,
        ccON_lower=on_result.cc_lower,
        ccON_median=on_result.cc_median,
        ccON_upper=on_result.cc_upper,
        ccON_se=on_result.cc_se,
        ccReporters_lower=reporter_result.cc_lower,
        ccReporters_median=reporter_result.cc_median,
        ccReporters_upper=reporter_result.cc_upper,
        ccReporters_se=reporter_result.cc_se,
        ac1_lower=isnothing(intensity_result) ? nothing : intensity_result.ac1_lower,
        ac1_median=isnothing(intensity_result) ? nothing : intensity_result.ac1_median,
        ac1_upper=isnothing(intensity_result) ? nothing : intensity_result.ac1_upper,
        ac1_se=isnothing(intensity_result) ? nothing : intensity_result.ac1_se,
        ac2_lower=isnothing(intensity_result) ? nothing : intensity_result.ac2_lower,
        ac2_median=isnothing(intensity_result) ? nothing : intensity_result.ac2_median,
        ac2_upper=isnothing(intensity_result) ? nothing : intensity_result.ac2_upper,
        ac2_se=isnothing(intensity_result) ? nothing : intensity_result.ac2_se,
        ac1ON_lower=on_result.ac1_lower,
        ac1ON_median=on_result.ac1_median,
        ac1ON_upper=on_result.ac1_upper,
        ac1ON_se=on_result.ac1_se,
        ac2ON_lower=on_result.ac2_lower,
        ac2ON_median=on_result.ac2_median,
        ac2ON_upper=on_result.ac2_upper,
        ac2ON_se=on_result.ac2_se,
        ac1Reporters_lower=reporter_result.ac1_lower,
        ac1Reporters_median=reporter_result.ac1_median,
        ac1Reporters_upper=reporter_result.ac1_upper,
        ac1Reporters_se=reporter_result.ac1_se,
        ac2Reporters_lower=reporter_result.ac2_lower,
        ac2Reporters_median=reporter_result.ac2_median,
        ac2Reporters_upper=reporter_result.ac2_upper,
        ac2Reporters_se=reporter_result.ac2_se
    )

    # Transform theory to match empirical based on CorrelationTrait
    # Theory: cc = E[XY] (raw uncentered correlation from covariance_functions)
    # Empirical: Uses CorrelationTrait directly - no transformations needed

    needs_centering = correlation_algorithm.centering != :none
    needs_normalization = correlation_algorithm.normalization != :none

    # Compute normalization factors (theoretical means)
    norm_cc = (m1 * m2) != 0 ? (m1 * m2) : 1.0
    norm_ac1 = m1^2 != 0 ? m1^2 : 1.0
    norm_ac2 = m2^2 != 0 ? m2^2 : 1.0
    norm_ccON = (m1ON * m2ON) != 0 ? (m1ON * m2ON) : 1.0
    norm_ac1ON = m1ON^2 != 0 ? m1ON^2 : 1.0
    norm_ac2ON = m2ON^2 != 0 ? m2ON^2 : 1.0
    norm_ccReporters = isnothing(m1Reporters) || isnothing(m2Reporters) || (m1Reporters * m2Reporters) != 0 ? (m1Reporters * m2Reporters) : 1.0
    norm_ac1Reporters = isnothing(m1Reporters) || m1Reporters^2 != 0 ? m1Reporters^2 : 1.0
    norm_ac2Reporters = isnothing(m2Reporters) || m2Reporters^2 != 0 ? m2Reporters^2 : 1.0

    # Step 1: Center (if centering is enabled)
    if needs_centering
        cc_centered = cc .- norm_cc
        ac1_centered = ac1 .- norm_ac1
        ac2_centered = ac2 .- norm_ac2
        ccON_centered = ccON .- norm_ccON
        ac1ON_centered = ac1ON .- norm_ac1ON
        ac2ON_centered = ac2ON .- norm_ac2ON
        ccReporters_centered = isnothing(ccReporters) ? nothing : ccReporters .- norm_ccReporters
        ac1Reporters_centered = isnothing(ac1Reporters) ? nothing : ac1Reporters .- norm_ac1Reporters
        ac2Reporters_centered = isnothing(ac2Reporters) ? nothing : ac2Reporters .- norm_ac2Reporters
    else
        cc_centered = cc
        ac1_centered = ac1
        ac2_centered = ac2
        ccON_centered = ccON
        ac1ON_centered = ac1ON
        ac2ON_centered = ac2ON
        ccReporters_centered = ccReporters
        ac1Reporters_centered = ac1Reporters
        ac2Reporters_centered = ac2Reporters
    end

    # Step 2: Normalize (if normalization is enabled)
    if needs_normalization
        cc_normalized = cc_centered ./ norm_cc
        ac1_normalized = ac1_centered ./ norm_ac1
        ac2_normalized = ac2_centered ./ norm_ac2
        ccON_normalized = ccON_centered ./ norm_ccON
        ac1ON_normalized = ac1ON_centered ./ norm_ac1ON
        ac2ON_normalized = ac2ON_centered ./ norm_ac2ON
        ccReporters_normalized = isnothing(ccReporters_centered) ? nothing : ccReporters_centered ./ norm_ccReporters
        ac1Reporters_normalized = isnothing(ac1Reporters_centered) ? nothing : ac1Reporters_centered ./ norm_ac1Reporters
        ac2Reporters_normalized = isnothing(ac2Reporters_centered) ? nothing : ac2Reporters_centered ./ norm_ac2Reporters
    else
        cc_normalized = cc_centered
        ac1_normalized = ac1_centered
        ac2_normalized = ac2_centered
        ccON_normalized = ccON_centered
        ac1ON_normalized = ac1ON_centered
        ac2ON_normalized = ac2ON_centered
        ccReporters_normalized = ccReporters_centered
        ac1Reporters_normalized = ac1Reporters_centered
        ac2Reporters_normalized = ac2Reporters_centered
    end

    # Theory autocovariances from `covariance_functions` are on positive lags only;
    # symmetrize to compare with empirical (which has symmetric lags from compute_covariance).
    # Note: cc, ccON, and ccReporters are already symmetrized in covariance_functions (lines 2717-2719),
    # but ac1, ac2, ac1ON, ac2ON, ac1Reporters, ac2Reporters are NOT symmetrized (they come from autocorfn_hmm/crosscorfn_hmm).
    # Use normalized versions for symmetrization
    ac1_theory_full = vcat(reverse(ac1_normalized), ac1_normalized[2:end])
    ac2_theory_full = vcat(reverse(ac2_normalized), ac2_normalized[2:end])
    ac1ON_theory_full = vcat(reverse(ac1ON_normalized), ac1ON_normalized[2:end])
    ac2ON_theory_full = vcat(reverse(ac2ON_normalized), ac2ON_normalized[2:end])
    ac1Reporters_theory_full = isnothing(ac1Reporters_normalized) ? nothing : vcat(reverse(ac1Reporters_normalized), ac1Reporters_normalized[2:end])
    ac2Reporters_theory_full = isnothing(ac2Reporters_normalized) ? nothing : vcat(reverse(ac2Reporters_normalized), ac2Reporters_normalized[2:end])
    # cc, ccON, and ccReporters are already symmetrized from covariance_functions, use normalized versions
    cc = cc_normalized
    ccON = ccON_normalized
    ccReporters = ccReporters_normalized

    # ============================================================================
    # FULL AUDIT: Validate theory vs empirical alignment
    # ============================================================================

    # Check lag alignment: theory uses full_lags (symmetric), empirical uses lags (should also be symmetric)
    n_lags_theory = length(lags)  # full_lags from covariance_functions
    n_lags_empirical = length(empirical.cc)

    # Verify that empirical lags match theory lags (they should be identical)
    # Check using on_result since it's always computed
    if hasproperty(on_result, :lags) && !isnothing(on_result.lags)
        empirical_lags = on_result.lags
        if length(empirical_lags) != length(lags) || !all(isapprox.(empirical_lags, lags, atol=1e-10))
            error("""
            LAG VALUES MISMATCH:
            Theory lags: $(lags[1:min(5, length(lags))]) ... $(lags[max(1, length(lags)-4):length(lags)])
            Empirical lags: $(empirical_lags[1:min(5, length(empirical_lags))]) ... $(empirical_lags[max(1, length(empirical_lags)-4):length(empirical_lags)])
            The lags passed to compute_covariance must exactly match full_lags from covariance_functions.
            """)
        end
    end

    if n_lags_theory != n_lags_empirical
        error("""
        LAG LENGTH MISMATCH:
        Theory has $(n_lags_theory) lags: $(lags[1:min(5, length(lags))]) ... $(lags[max(1, length(lags)-4):length(lags)])
        Empirical has $(n_lags_empirical) lags
        Check that lags passed to compute_correlation_functions match full_lags from correlation_functions.
        """)
    end

    # Validate that lags are symmetric (theory should have symmetric lags)
    if length(lags) > 1
        mid_idx = div(length(lags), 2) + 1
        if lags[mid_idx] != 0.0
            @warn "Zero lag not at center: lags[$mid_idx] = $(lags[mid_idx]), expected 0.0"
        end
        # Check symmetry: lags should be symmetric around zero
        if length(lags) > 1
            neg_lags = lags[1:(mid_idx-1)]
            pos_lags = lags[(mid_idx+1):end]
            if length(neg_lags) == length(pos_lags)
                if !all(isapprox.(-reverse(pos_lags), neg_lags, atol=1e-10))
                    @warn "Lags are not symmetric around zero"
                end
            end
        end
    end

    # Validate all theory quantities have correct lengths
    @assert length(cc) == n_lags_theory "cc theory length mismatch: $(length(cc)) != $n_lags_theory"
    @assert length(ccON) == n_lags_theory "ccON theory length mismatch: $(length(ccON)) != $n_lags_theory"
    @assert length(ac1_theory_full) == n_lags_theory "ac1_theory_full length mismatch: $(length(ac1_theory_full)) != $n_lags_theory"
    @assert length(ac2_theory_full) == n_lags_theory "ac2_theory_full length mismatch: $(length(ac2_theory_full)) != $n_lags_theory"
    @assert length(ac1ON_theory_full) == n_lags_theory "ac1ON_theory_full length mismatch: $(length(ac1ON_theory_full)) != $n_lags_theory"
    @assert length(ac2ON_theory_full) == n_lags_theory "ac2ON_theory_full length mismatch: $(length(ac2ON_theory_full)) != $n_lags_theory"

    # Validate all empirical quantities have correct lengths
    @assert length(empirical.cc) == n_lags_empirical "empirical.cc length mismatch"
    @assert length(empirical.ac1) == n_lags_empirical "empirical.ac1 length mismatch"
    @assert length(empirical.ac2) == n_lags_empirical "empirical.ac2 length mismatch"
    @assert length(empirical.ccON) == n_lags_empirical "empirical.ccON length mismatch"
    @assert length(empirical.ac1ON) == n_lags_empirical "empirical.ac1ON length mismatch"
    @assert length(empirical.ac2ON) == n_lags_empirical "empirical.ac2ON length mismatch"
    @assert length(empirical.ccReporters) == n_lags_empirical "empirical.ccReporters length mismatch"
    @assert length(empirical.ac1Reporters) == n_lags_empirical "empirical.ac1Reporters length mismatch"
    @assert length(empirical.ac2Reporters) == n_lags_empirical "empirical.ac2Reporters length mismatch"

    # Validate reporter theory lengths if provided
    if !isnothing(ccReporters)
        @assert length(ccReporters) == n_lags_theory "ccReporters theory length mismatch: $(length(ccReporters)) != $n_lags_theory"
        @assert length(ac1Reporters_theory_full) == n_lags_theory "ac1Reporters_theory_full length mismatch: $(length(ac1Reporters_theory_full)) != $n_lags_theory"
        @assert length(ac2Reporters_theory_full) == n_lags_theory "ac2Reporters_theory_full length mismatch: $(length(ac2Reporters_theory_full)) != $n_lags_theory"
    end

    # Check autocovariance symmetry: ac1, ac2, ac1ON, ac2ON should be symmetric (ac(τ) = ac(-τ))
    # Theory autocovariances after symmetrization
    mid_idx = div(n_lags_theory, 2) + 1
    @assert ac1_theory_full[mid_idx] ≈ ac1_theory_full[mid_idx] "ac1 should have same value at zero lag"
    if n_lags_theory > 3
        # Check that ac1 is symmetric: ac1[-τ] should equal ac1[τ]
        neg_side = ac1_theory_full[1:(mid_idx-1)]
        pos_side = ac1_theory_full[(mid_idx+1):end]
        if !all(isapprox.(reverse(neg_side), pos_side, atol=1e-10))
            @warn "ac1_theory_full is not symmetric around zero lag"
        end
    end

    # Compute norms comparing theory and empirical (both should have same length: symmetric lags)
    linf_norm = [maximum(abs.(cc .- empirical.cc))]
    l2_norm = [sqrt(sum((cc .- empirical.cc) .^ 2))]

    # Diagnostic: Print summary statistics for theory vs empirical comparison
    # Check zero-lag values for ON state autocovariance (should match if uncentered)
    # Theory returns uncentered R_ON_ON(0) = E[ON²], empirical returns uncentered R_ON_ON(0) = E[ON²]
    # At zero lag, both should equal the second moment E[ON²] = mON (since ON is binary, ON² = ON)
    @info """
    THEORY VS EMPIRICAL AUDIT SUMMARY:
    =================================
    Lag range: [$(lags[1]), ..., $(lags[end])] ($(n_lags_theory) lags)
    Number of trials: $(ntrials)
    Trace lengths: $(minimum([size(t, 1) for t in intensity_traces])) - $(maximum([size(t, 1) for t in intensity_traces])) frames

    ON state means (theory vs empirical):
      Theory mON1: $(m1ON), mON2: $(m2ON)
      Empirical mON1: $(empirical.mON1), mON2: $(empirical.mON2)
      Mean difference unit 1: $(m1ON - empirical.mON1) ($((m1ON - empirical.mON1) / m1ON * 100)% relative)
      Mean difference unit 2: $(m2ON - empirical.mON2) ($((m2ON - empirical.mON2) / m2ON * 100)% relative)
      Expected zero-lag for ac1ON (theory): $(m1ON) (since ON is binary, E[ON²] = E[ON] = mON)
      Expected zero-lag for ac1ON (empirical): $(empirical.mON1) (since ON is binary, E[ON²] = E[ON])
      Expected zero-lag for ac2ON (theory): $(m2ON) (since ON is binary, E[ON²] = E[ON] = mON)
      Expected zero-lag for ac2ON (empirical): $(empirical.mON2) (since ON is binary, E[ON²] = E[ON])

    Intensity cross-covariance (cc):
      Theory zero-lag: $(cc[mid_idx])
      Empirical zero-lag (mean): $(empirical.cc[mid_idx])
      L∞ norm: $(linf_norm[1])
      L² norm: $(l2_norm[1])

    ON state cross-covariance (ccON):
      Theory zero-lag: $(ccON[mid_idx])
      Empirical zero-lag (mean): $(empirical.ccON[mid_idx])
      Theory max: $(maximum(ccON)) at lag $(lags[argmax(ccON)])
      Empirical max: $(maximum(empirical.ccON)) at lag $(lags[argmax(empirical.ccON)])

    ON state autocovariance unit 1 (ac1ON):
      Theory zero-lag: $(ac1ON_theory_full[mid_idx])
      Empirical zero-lag: $(empirical.ac1ON[mid_idx])
      Difference at zero-lag: $(ac1ON_theory_full[mid_idx] - empirical.ac1ON[mid_idx])
      Relative difference at zero-lag: $((ac1ON_theory_full[mid_idx] - empirical.ac1ON[mid_idx]) / ac1ON_theory_full[mid_idx] * 100)%
      Theory max: $(maximum(ac1ON_theory_full)) at lag $(lags[argmax(ac1ON_theory_full)])
      Empirical max: $(maximum(empirical.ac1ON)) at lag $(lags[argmax(empirical.ac1ON)])

    ON state autocovariance unit 2 (ac2ON):
      Theory zero-lag: $(ac2ON_theory_full[mid_idx])
      Empirical zero-lag (mean): $(empirical.ac2ON[mid_idx])
      Difference at zero-lag: $(ac2ON_theory_full[mid_idx] - empirical.ac2ON[mid_idx])
      Relative difference at zero-lag: $((ac2ON_theory_full[mid_idx] - empirical.ac2ON[mid_idx]) / ac2ON_theory_full[mid_idx] * 100)%
      Theory max: $(maximum(ac2ON_theory_full)) at lag $(lags[argmax(ac2ON_theory_full)])
      Empirical max: $(maximum(empirical.ac2ON)) at lag $(lags[argmax(empirical.ac2ON)])
    """

    return (
        # Intensity correlations (unnormalized)
        cc_theory=cc, cc_mean=empirical.cc,
        cc_se=empirical.cc_se,
        cc_lower=empirical.cc_lower,
        cc_median=empirical.cc_median,
        cc_upper=empirical.cc_upper,
        ac1_mean=empirical.ac1, ac1_theory=ac1_theory_full,
        ac1_lower=empirical.ac1_lower, ac1_median=empirical.ac1_median, ac1_upper=empirical.ac1_upper, ac1_se=empirical.ac1_se,
        ac2_mean=empirical.ac2, ac2_theory=ac2_theory_full,
        ac2_lower=empirical.ac2_lower, ac2_median=empirical.ac2_median, ac2_upper=empirical.ac2_upper, ac2_se=empirical.ac2_se,
        # ON state correlations (unnormalized)
        ccON_theory=ccON, ccON_mean=empirical.ccON,
        ccON_se=empirical.ccON_se,
        ccON_lower=empirical.ccON_lower,
        ccON_median=empirical.ccON_median,
        ccON_upper=empirical.ccON_upper,
        ac1ON_mean=empirical.ac1ON, ac1ON_theory=ac1ON_theory_full,
        ac1ON_lower=empirical.ac1ON_lower, ac1ON_median=empirical.ac1ON_median, ac1ON_upper=empirical.ac1ON_upper, ac1ON_se=empirical.ac1ON_se,
        ac2ON_mean=empirical.ac2ON, ac2ON_theory=ac2ON_theory_full,
        ac2ON_lower=empirical.ac2ON_lower, ac2ON_median=empirical.ac2ON_median, ac2ON_upper=empirical.ac2ON_upper, ac2ON_se=empirical.ac2ON_se,
        # Reporter count correlations (unnormalized)
        ccReporters_theory=ccReporters, ccReporters_mean=empirical.ccReporters,
        ccReporters_se=empirical.ccReporters_se,
        ccReporters_lower=empirical.ccReporters_lower,
        ccReporters_median=empirical.ccReporters_median,
        ccReporters_upper=empirical.ccReporters_upper,
        ac1Reporters_mean=empirical.ac1Reporters, ac1Reporters_theory=ac1Reporters_theory_full,
        ac1Reporters_lower=empirical.ac1Reporters_lower, ac1Reporters_median=empirical.ac1Reporters_median, ac1Reporters_upper=empirical.ac1Reporters_upper, ac1Reporters_se=empirical.ac1Reporters_se,
        ac2Reporters_mean=empirical.ac2Reporters, ac2Reporters_theory=ac2Reporters_theory_full,
        ac2Reporters_lower=empirical.ac2Reporters_lower, ac2Reporters_median=empirical.ac2Reporters_median, ac2Reporters_upper=empirical.ac2Reporters_upper, ac2Reporters_se=empirical.ac2Reporters_se,
        # Other statistics
        linf_norm=linf_norm, l2_norm=l2_norm,
        lags=lags,  # Symmetric lags (full_lags from covariance_functions) - same for theory and empirical
        mean1=empirical.mON1, mean2=empirical.mON2,
        m1=m1, m2=m2, v1=v1, v2=v2, v1_empirical=empirical.v1_empirical, v2_empirical=empirical.v2_empirical,
        m1ON=m1ON, m2ON=m2ON, mR1_empirical=empirical.mR1, mR2_empirical=empirical.mR2, m1Reporters=m1Reporters, m2Reporters=m2Reporters,
        # Diagnostic info for validation
        n_lags=n_lags_theory,  # Number of lags (should be same for all theory and empirical)
        trace_lengths=[size(t, 1) for t in intensity_traces],  # Length of each trace (for validation)
        n_trials=ntrials  # Number of simulation trials
    )
end


"""
    write_correlation_functions_file(file, transitions=..., G=(3, 3), R=(3, 3), S=(1, 0), insertstep=(1, 1), pattern="gene", lags=collect(0:10:1000), probfn=prob_Gaussian, ratetype="ml")

Compute theoretical correlation functions for a single rate file and return the results.

# Arguments
- `file::String`: Path to the rate file (typically a `.txt` file containing rate parameters)
- `transitions::Tuple`: Tuple of transition definitions for each unit (default: 3-state model transitions)
- `G::Tuple`: Number of gene states for each unit (default: (3, 3))
- `R::Tuple`: Number of RNA states for each unit (default: (3, 3))
- `S::Tuple`: Initial state definitions (default: (1, 0))
- `insertstep::Tuple`: Insert step definitions (default: (1, 1))
- `interval::Float64`: Time interval for transitions (default: 1.0)
- `pattern::String`: Pattern to extract coupling information from filename (default: "gene")
- `lags::Vector{Int}`: Time lags for correlation function computation (default: 0:10:1000)
- `probfn`: Probability function for observation model (default: prob_Gaussian)
- `ratetype::String`: Which row of rates to use from file (default: "ml" for maximum likelihood)

# Returns
- `Tuple` containing:
  - `ac1, ac2`: Auto-correlation functions for intensity (unit1, unit2)
  - `cc`: Cross-correlation function for intensity
  - `ccON`: Cross-correlation function for ON states (unnormalized: E[xy] - E[x]E[y])
  - `tau`: Time lags (with negative lags included)
  - `m1, m2`: Mean intensities
  - `v1, v2`: Variances
  - `mON1, mON2`: Mean ON state probabilities
  - `ac1ON, ac2ON`: Auto-correlation functions for ON states
  - `ccReporters`: Cross-correlation function for reporter counts (unnormalized)
  - `mR1, mR2`: Mean reporter counts
  - `ac1Reporters, ac2Reporters`: Auto-correlation functions for reporter counts

# Notes
- Extracts coupling information from the filename using `pattern`
- All cross-correlations and auto-correlations are unnormalized (E[xy] - E[x]E[y])
- ON states are binary (1 if reporter count > 0, 0 otherwise)
"""
function write_correlation_functions_file(file, transitions=(([1, 2], [2, 1], [2, 3], [3, 2]), ([1, 2], [2, 1], [2, 3], [3, 2])), G=(3, 3), R=(3, 3), S=(1, 0), insertstep=(1, 1), pattern="gene", lags=collect(0:1:200), probfn=prob_Gaussian, ratetype="ml")
    println(file)
    r = readrates(file, get_row(ratetype))
    source, target = extract_source_target(pattern, file)
    (source == "R") && (source = collect(G[1]+1:G[1]+R[1]))
    coupling = ((1, 2), (tuple(), tuple(1)), (source, 0), (0, target), 1)
    correlation_functions(r, transitions, G, R, S, insertstep, probfn, coupling, lags)
end

"""
    write_correlation_functions(folder, transitions=..., G=(3, 3), R=(3, 3), S=(0, 0), insertstep=(1, 1), pattern="gene", lags=collect(0:1:200), probfn=prob_Gaussian, ratetype="median")

Compute theoretical correlation functions for all rate files in a folder and write results to CSV files.

This is the main function for batch-processing rate parameter files to generate theoretical cross-correlation
and auto-correlation predictions. It recursively searches a folder for files matching the pattern
`*rates*tracejoint*`, computes theoretical correlation functions for each using the HMM framework,
and writes comprehensive results to CSV files that can be used for model scoring against empirical data.

# Arguments
- `folder::String`: Path to folder containing rate files. The function recursively searches all subdirectories.
- `transitions::Tuple`: Tuple of transition definitions for each unit. Each element specifies allowed
  state transitions as a vector of pairs `[from, to]`. Default: `(([1, 2], [2, 1], [2, 3], [3, 2]), ([1, 2], [2, 1], [2, 3], [3, 2]))`
  for a 3-state model with forward/backward transitions for both units.
- `G::Tuple{Int, Int}`: Number of gene states for each unit (default: `(3, 3)`).
- `R::Tuple{Int, Int}`: Number of RNA states for each unit (default: `(3, 3)`).
- `S::Tuple{Int, Int}`: Initial state definitions (default: `(0, 0)`).
- `insertstep::Tuple{Int, Int}`: Insert step definitions (default: `(1, 1)`).
- `interval::Float64`: Time interval for transitions in the same units as the rate parameters (default: `1.0`).
- `pattern::String`: Pattern to extract coupling information from filename (default: `"gene"`). Used by
  `extract_source_target` to determine which states are coupled.
- `lags::Vector{Int}`: Time lags for correlation function computation. The function computes correlations at these
  lags and also includes negative lags (symmetric). Default: `collect(0:1:200)` which gives lags from -200 to 200.
- `probfn`: Probability function for the observation model (default: `prob_Gaussian`). Determines how
  reporter counts are generated from hidden states.
- `ratetype::String`: Which row of rates to use from the rate file (default: `"median"`). Other options
  include `"ml"` (maximum likelihood), `"mean"`, etc. This selects which parameter estimate to use when
  multiple estimates are stored in the file.

# Output Files

For each input file matching `*rates*tracejoint*.txt`, creates a corresponding output file
`*crosscorrelation*tracejoint*.csv` in the same directory. The CSV file contains the following columns:

## Time Lags
- `tau::Vector{Int}`: Time lags from `-max_lag` to `+max_lag` (symmetric around zero)

## ON State Correlation Functions (Binary: 1 if reporter > 0, 0 otherwise)
- `cc_ON::Vector{Float64}`: Cross-correlation function between enhancer and gene ON states (unnormalized: E[xy] - E[x]E[y]).
  Positive τ means enhancer leads (E[enhancer(t) × gene(t+τ)] - E[enhancer] × E[gene]).
- `ac1_ON::Vector{Float64}`: Auto-correlation function of enhancer ON states (unnormalized, symmetric: includes negative lags).
- `ac2_ON::Vector{Float64}`: Auto-correlation function of gene ON states (unnormalized, symmetric: includes negative lags).
- `mON1::Vector{Float64}`: Mean ON state probability for enhancer (repeated for each lag, scalar value).
- `mON2::Vector{Float64}`: Mean ON state probability for gene (repeated for each lag, scalar value).

## Reporter Count Correlation Functions (Raw integer counts)
- `cc_Reporters::Vector{Float64}`: Cross-correlation function between enhancer and gene reporter counts (unnormalized).
  Same convention as `cc_ON` (positive τ means enhancer leads).
- `ac1_Reporters::Vector{Float64}`: Auto-correlation function of enhancer reporter counts (unnormalized, symmetric).
- `ac2_Reporters::Vector{Float64}`: Auto-correlation function of gene reporter counts (unnormalized, symmetric).
- `mR1::Vector{Float64}`: Mean reporter count for enhancer (repeated for each lag, scalar value).
- `mR2::Vector{Float64}`: Mean reporter count for gene (repeated for each lag, scalar value).

# Workflow

1. **File Discovery**: Recursively walks through `folder` and identifies all files containing both
   `"rates"` and `"tracejoint"` in the filename.

2. **Rate Loading**: For each matching file, loads rate parameters using `readrates(file, get_row(ratetype))`.

3. **Coupling Extraction**: Extracts coupling information from the filename using `extract_source_target(pattern, file)`
   to determine which states are coupled (e.g., gene state 1 to gene state 3, or RNA states).

4. **Correlation Function Computation**: Calls `correlation_functions` (via `write_correlation_functions_file`) to compute all theoretical
   correlation functions using the HMM framework.

5. **File Writing**: Writes results to CSV with filename pattern: `rates_*.txt` → `crosscorrelation_*.csv`.

# Usage Example

```julia
# Process all rate files in a results folder
write_correlation_functions(
    "results/5Prime-coupled-2025-11-27/",
    transitions=(([1, 2], [2, 1], [2, 3], [3, 2]), ([1, 2], [2, 1], [2, 3], [3, 2])),
    G=(3, 3),
    R=(3, 3),
    lags=collect(0:1:200),  # Lags from -200 to 200
    ratetype="median"
)
```

# Notes

- **Unnormalized Correlation Functions**: All correlation function values are unnormalized (E[xy] - E[x]E[y]). To normalize
  by means, divide by (m1 × m2) for cross-correlations or by (m × m) for auto-correlations.

- **Lag Convention**: Positive τ means the first unit (enhancer) leads the second unit (gene). This matches
  the convention used in `StatsBase.crosscov(enhancer, gene, lags)`.

- **Symmetric Auto-correlations**: Auto-correlation functions are symmetric (ac(τ) = ac(-τ)), so negative lags are
  included by reversing the positive lag values.

- **File Matching**: Only files containing both `"rates"` and `"tracejoint"` in the filename are processed.
  This pattern identifies files containing joint trace fitting results.

- **Batch Processing**: This function is designed for batch processing of multiple model fits. Each rate
  file typically corresponds to a different coupling model (e.g., gene state 1→3 coupling vs. gene state 2→3).

- **Model Scoring**: The output CSV files are designed to be read by `score_models_from_traces`, which
  compares these theoretical predictions against empirical correlation functions computed from experimental
  or simulated trace data.

# See Also
- `write_correlation_functions_file`: Processes a single rate file (called internally by this function)
- `correlation_functions`: Core HMM function that computes theoretical correlation functions
- `score_models_from_traces`: Scores theoretical predictions against empirical data
- `readrates`: Loads rate parameters from text files
"""
function write_correlation_functions(folder; transitions=(([1, 2], [2, 1], [2, 3], [3, 2]), ([1, 2], [2, 1], [2, 3], [3, 2])), G=(3, 3), R=(3, 3), S=(0, 0), insertstep=(1, 1), pattern="gene", lags=collect(0:1:200), probfn=prob_Gaussian, ratetype="median")
    for (root, dirs, files) in walkdir(folder)
        for f in files
            if occursin("rates", f) && occursin("tracejoint", f)
                file = joinpath(root, f)
                tau, cc, ac1, ac2, m1, m2, v1, v2, ccON, ac1ON, ac2ON, mON1, mON2, v1ON, v2ON, ccReporters, ac1Reporters, ac2Reporters, mReporters1, mReporters2, v1Reporters, v2Reporters = write_correlation_functions_file(file, transitions, G, R, S, insertstep, pattern, lags, probfn, ratetype)
                parts = fields(basename(file))
                new_model = create_modelstring(G, R, S, insertstep)
                out = joinpath(root, "crosscorrelation_" * parts.label * "_" * parts.cond * "_" * parts.gene * "_" * new_model * "_" * parts.nalleles * ".csv")
                n_lags = length(tau)

                CSV.write(out, DataFrame(
                    tau=tau,
                    # ON state: unnormalized cc, ac1, ac2, m1, m2
                    cc_ON=ccON,
                    ac1_ON=[reverse(ac1ON); ac1ON[2:end]],
                    ac2_ON=[reverse(ac2ON); ac2ON[2:end]],
                    m_ON1=fill(mON1, n_lags),
                    m_ON2=fill(mON2, n_lags),
                    # Reporters: unnormalized cc, ac1, ac2, m1, m2
                    cc_Reporters=ccReporters,
                    ac1_Reporters=[reverse(ac1Reporters); ac1Reporters[2:end]],
                    ac2_Reporters=[reverse(ac2Reporters); ac2Reporters[2:end]],
                    m_Reporters1=fill(mReporters1, n_lags),
                    mReporters2=fill(mReporters2, n_lags)
                ))
            end
        end
    end
end

"""
    write_correlation_functions_empirical(trace_folder, unit1_filename, unit2_filename, output_file; intensity_labels=nothing, intensity_start=1, intensity_stop=-1, intensity_col=3, lags=collect(0:1:200), on_threshold=0, mON1_theory=nothing, mON2_theory=nothing, mR1_theory=nothing, mR2_theory=nothing, bootstrap=false, n_bootstrap=1000)

Compute empirical correlation functions from trace files and write results to CSV file.

This function computes cross-correlation functions for:
1. **Intensity**: From experimental intensity traces (if `intensity_datapath` and `intensity_labels` provided)
2. **Predicted ON states**: From predicted trace files (binary: 1 if reporter > threshold, 0 otherwise)
3. **Predicted reporter counts**: From predicted trace files (raw reporter counts)

Then writes the results to a CSV file in the same format as `write_correlation_functions` output.

# Arguments
- `trace_folder::String`: Path to folder containing intensity trace files
- `unit1_filename::String`: Path to predicted trace file for unit 1 (enhancer) (from uncoupled model)
- `unit2_filename::String`: Path to predicted trace file for unit 2 (gene) (from uncoupled model)
- `output_file::String`: Path to output CSV file
- `intensity_labels::Union{Vector{String}, Nothing}=nothing`: Labels to match in filenames for joint traces (e.g., `["enhancer", "gene"]`). The function will find files containing each label in their filename and pair them (typically by numeric prefix, e.g., `001_enhancer.trk` pairs with `001_gene.trk`). Each pair becomes a joint trace matrix with columns `[enhancer, gene]`. If `nothing`, intensity traces are not loaded.
- `traceinfo::Tuple=(1.0, 1.0, -1)`: Trace parameters tuple `(dt, start, stop)` where `dt` is the frame interval in minutes, `start` is the starting frame time in minutes, and `stop` is the ending frame time (-1 for last frame). Same format as used in `load_data_trace`.
- `intensity_col::Int=3`: Column to read from intensity trace files (default: 3)
- `lags::Vector{Int}=collect(0:1:200)`: Time lags for cross-correlation computation
- `on_threshold::Real=0`: Threshold for converting reporter counts to ON/OFF states
- `mON1_theory::Union{Float64, Nothing}=nothing`: Theoretical mean ON state for enhancer (if provided, used for centered version)
- `mON2_theory::Union{Float64, Nothing}=nothing`: Theoretical mean ON state for gene (if provided, used for centered version)
- `mR1_theory::Union{Float64, Nothing}=nothing`: Theoretical mean reporter count for enhancer (if provided, used for centered version)
- `mR2_theory::Union{Float64, Nothing}=nothing`: Theoretical mean reporter count for gene (if provided, used for centered version)
- `bootstrap::Bool=false`: Whether to compute bootstrap confidence intervals
- `n_bootstrap::Int=1000`: Number of bootstrap replicates (if bootstrap=true)
- `correlation_algorithm::CorrelationTrait=StandardCorrelation()`: Correlation algorithm to use. Options include `StandardCorrelation()`, `WindowedCorrelation()`, `MultiTauCorrelation()`, `IDLCorrelation()`, or custom `CorrelationTrait(centering=:global_mean, ...)`

# Output CSV Structure

The CSV file contains the same structure as `write_correlation_functions` output:

## Time Lags
- `tau::Vector{Int}`: Time lags from `-max_lag` to `+max_lag` (symmetric around zero)

## ON State Correlation Functions (Binary: 1 if reporter > threshold, 0 otherwise)
- `cc_ON::Vector{Float64}`: Empirical cross-correlation function between enhancer and gene ON states (unnormalized).
- `ac1_ON::Vector{Float64}`: Empirical auto-correlation function of enhancer ON states (unnormalized, symmetric).
- `ac2_ON::Vector{Float64}`: Empirical auto-correlation function of gene ON states (unnormalized, symmetric).
- `mON1::Vector{Float64}`: Mean ON state probability for enhancer (repeated for each lag).
- `mON2::Vector{Float64}`: Mean ON state probability for gene (repeated for each lag).

## Reporter Count Correlation Functions (Raw integer counts)
- `cc_Reporters::Vector{Float64}`: Empirical cross-correlation function between enhancer and gene reporter counts (unnormalized).
- `ac1_Reporters::Vector{Float64}`: Empirical auto-correlation function of enhancer reporter counts (unnormalized, symmetric).
- `ac2_Reporters::Vector{Float64}`: Empirical auto-correlation function of gene reporter counts (unnormalized, symmetric).
- `mR1::Vector{Float64}`: Mean reporter count for enhancer (repeated for each lag).
- `mR2::Vector{Float64}`: Mean reporter count for gene (repeated for each lag).

## Intensity Correlation Functions (if intensity traces provided)
- `cc::Vector{Float64}`: Empirical cross-correlation function between enhancer and gene intensity (unnormalized).
- `ac1::Vector{Float64}`: Empirical auto-correlation function of enhancer intensity (unnormalized, symmetric).
- `ac2::Vector{Float64}`: Empirical auto-correlation function of gene intensity (unnormalized, symmetric).
- `m1::Vector{Float64}`: Mean intensity for enhancer (repeated for each lag).
- `m2::Vector{Float64}`: Mean intensity for gene (repeated for each lag).

If `bootstrap=true`, also includes confidence intervals (`_lower`, `_median`, `_upper`, `_se`) for all quantities.

# Algorithm

1. **Compute Cross-Correlations**: Calls `compute_correlation_functions` to compute empirical cross-correlation functions
2. **Symmetrize Lags**: Adds negative lags for symmetric output (same as `write_correlation_functions`)
3. **Write CSV**: Writes results to CSV file in the same format as `write_correlation_functions`

# Notes

- The output format matches `write_correlation_functions` exactly, ensuring direct comparability.
- Uses the same lags as specified, with symmetric negative lags included.
- All correlation functions are unnormalized (centered: E[xy] - E[x]E[y]).

# Usage Example

```julia
# Compute empirical correlation functions from predicted trace files only (no intensity)
write_correlation_functions_empirical(
    "data/traces/",
    "results/predicted_enhancer.csv",
    "results/predicted_gene.csv",
    "results/crosscovariance_empirical.csv",
    lags=collect(0:1:200),
    bootstrap=true,
    n_bootstrap=1000
)

# With intensity traces from folder (files like 001_enhancer.trk, 001_gene.trk, etc.)
write_correlation_functions_empirical(
    "data/3Prime_gene_enhancer/including_background/short/",
    "results/predicted_enhancer.csv",
    "results/predicted_gene.csv",
    "results/crosscovariance_empirical.csv",
    intensity_labels=["enhancer", "gene"],  # Matches files: 001_enhancer.trk with 001_gene.trk, etc.
    traceinfo=(1.0, 1.0, -1),  # (dt, start, stop) in minutes, -1 means last frame
    intensity_col=3,
    lags=collect(0:1:200),
    bootstrap=true
)
```
"""
function write_correlation_functions_empirical(trace_folder::String, unit1_filename::String, unit2_filename::String, output_file::String; intensity_labels=nothing, traceinfo=(1.0, 1.0, -1), intensity_col=3, lags=collect(0:1:200), on_threshold=0, mON1_theory=nothing, mON2_theory=nothing, mR1_theory=nothing, mR2_theory=nothing, bootstrap=false, n_bootstrap=1000, correlation_algorithm=StandardCorrelation())
    # Load intensity traces if labels provided (following load_data_trace pattern)
    intensity_traces = Vector{Matrix{Float64}}()
    if !isempty(trace_folder)
        if !isnothing(intensity_labels)
            if intensity_labels isa Vector{String}
                # Joint traces: use read_tracefiles with traceinfo tuple format (like load_data_trace)
                traces_raw = read_tracefiles(trace_folder, intensity_labels, traceinfo, intensity_col)
                if length(traces_raw) == 0
                    throw("No traces found for intensity_labels=$intensity_labels in $trace_folder")
                end
                # Apply zero_median with zeromedian=true (following load_data_trace pattern)
                traces_zeroed, _ = zero_median(traces_raw, true)
                # Convert to Vector{Matrix{Float64}} to match expected type
                intensity_traces = Vector{Matrix{Float64}}(map(t -> Matrix{Float64}(Float64.(t)), traces_zeroed))
            else
                error("intensity_labels must be a Vector{String} for joint traces (e.g., [\"enhancer\", \"gene\"])")
            end
        end
    end

    # Check if we have intensity results
    has_intensity = !isnothing(intensity_labels) && !isempty(intensity_traces)

    # Load predicted traces (from uncoupled models) - gives reporter counts and ON states
    predicted_traces = prepare_traces_from_prediction(unit1_filename, unit2_filename; on_threshold=on_threshold)
    reporter_traces = predicted_traces.reporter_traces
    on_traces = predicted_traces.on_traces

    # Construct symmetric lags (same as write_correlation_functions): tau = [-max_lag, ..., -1, 0, 1, ..., max_lag]
    # Pass symmetric lags to compute_correlation_functions - correlation_function handles negative lags correctly
    lags_symmetric = vcat(-reverse(lags[2:end]), lags)
    tau = lags_symmetric
    n_lags = length(tau)

    # Compute empirical cross-correlation functions with symmetric lags
    # Use method that computes means from traces (most common case)
    if !isnothing(mON1_theory) && !isnothing(mON2_theory)
        on_result_full = compute_correlation_functions(on_traces, lags_symmetric; mean1=mON1_theory, mean2=mON2_theory, bootstrap=bootstrap, n_bootstrap=n_bootstrap, correlation_algorithm=correlation_algorithm)
    else
        on_result_full = compute_correlation_functions(on_traces, lags_symmetric; bootstrap=bootstrap, n_bootstrap=n_bootstrap, correlation_algorithm=correlation_algorithm)
    end
    if !isnothing(mR1_theory) && !isnothing(mR2_theory)
        reporter_result_full = compute_correlation_functions(reporter_traces, lags_symmetric; mean1=mR1_theory, mean2=mR2_theory, bootstrap=bootstrap, n_bootstrap=n_bootstrap, correlation_algorithm=correlation_algorithm)
    else
        reporter_result_full = compute_correlation_functions(reporter_traces, lags_symmetric; bootstrap=bootstrap, n_bootstrap=n_bootstrap, correlation_algorithm=correlation_algorithm)
    end

    # Results already have both positive and negative lags from crosscorrelation_function
    cc_ON_full = on_result_full.cc
    cc_Reporters_full = reporter_result_full.cc

    # Autocovariances ARE symmetric, but crosscorrelation_function already handles negative lags correctly
    ac1_ON_sym = on_result_full.ac1
    ac2_ON_sym = on_result_full.ac2
    ac1_Reporters_sym = reporter_result_full.ac1
    ac2_Reporters_sym = reporter_result_full.ac2

    # Build output DataFrame with same column order as write_correlation_functions
    # write_correlation_functions order: tau, cc_ON, ac1_ON, ac2_ON, m_ON1, m_ON2, cc_Reporters, ac1_Reporters, ac2_Reporters, m_Reporters1, mReporters2
    # Use OrderedDict or construct DataFrame with columns in correct order
    df = DataFrame(
        tau=tau,
        cc_ON=cc_ON_full,
        ac1_ON=ac1_ON_sym,
        ac2_ON=ac2_ON_sym,
        m_ON1=fill(on_result_full.mean1, n_lags),
        m_ON2=fill(on_result_full.mean2, n_lags),
        cc_Reporters=cc_Reporters_full,
        ac1_Reporters=ac1_Reporters_sym,
        ac2_Reporters=ac2_Reporters_sym,
        m_Reporters1=fill(reporter_result_full.mean1, n_lags),
        mReporters2=fill(reporter_result_full.mean2, n_lags)
    )

    # Add intensity correlation functions if available (after base columns, same order as write_correlation_functions would if it had intensity)
    if has_intensity
        # Compute intensity with symmetric lags as well (already includes both positive and negative lags)
        intensity_result_full = compute_correlation_functions(intensity_traces, lags_symmetric; bootstrap=bootstrap, n_bootstrap=n_bootstrap, correlation_algorithm=correlation_algorithm)
        cc_intensity_full = intensity_result_full.cc
        ac1_intensity_sym = intensity_result_full.ac1  # Already includes both positive and negative lags
        ac2_intensity_sym = intensity_result_full.ac2  # Already includes both positive and negative lags
        df.cc = cc_intensity_full
        df.ac1 = ac1_intensity_sym
        df.ac2 = ac2_intensity_sym
        df.m1 = fill(intensity_result_full.mean1, n_lags)
        df.m2 = fill(intensity_result_full.mean2, n_lags)
    end

    # Add bootstrap fields if available (after all base columns)
    # Bootstrap results already include both positive and negative lags
    if bootstrap && hasproperty(on_result_full, :cc_lower)
        # ON state bootstrap fields
        df.cc_ON_lower = on_result_full.cc_lower
        df.cc_ON_median = on_result_full.cc_median
        df.cc_ON_upper = on_result_full.cc_upper
        df.cc_ON_se = on_result_full.cc_se
        df.ac1_ON_lower = on_result_full.ac1_lower
        df.ac1_ON_median = on_result_full.ac1_median
        df.ac1_ON_upper = on_result_full.ac1_upper
        df.ac1_ON_se = on_result_full.ac1_se
        df.ac2_ON_lower = on_result_full.ac2_lower
        df.ac2_ON_median = on_result_full.ac2_median
        df.ac2_ON_upper = on_result_full.ac2_upper
        df.ac2_ON_se = on_result_full.ac2_se

        # Reporter bootstrap fields
        df.cc_Reporters_lower = reporter_result_full.cc_lower
        df.cc_Reporters_median = reporter_result_full.cc_median
        df.cc_Reporters_upper = reporter_result_full.cc_upper
        df.cc_Reporters_se = reporter_result_full.cc_se
        df.ac1_Reporters_lower = reporter_result_full.ac1_lower
        df.ac1_Reporters_median = reporter_result_full.ac1_median
        df.ac1_Reporters_upper = reporter_result_full.ac1_upper
        df.ac1_Reporters_se = reporter_result_full.ac1_se
        df.ac2_Reporters_lower = reporter_result_full.ac2_lower
        df.ac2_Reporters_median = reporter_result_full.ac2_median
        df.ac2_Reporters_upper = reporter_result_full.ac2_upper
        df.ac2_Reporters_se = reporter_result_full.ac2_se

        # Add intensity bootstrap fields if available
        if has_intensity && hasproperty(intensity_result_full, :cc_lower)
            df.cc_lower = intensity_result_full.cc_lower
            df.cc_median = intensity_result_full.cc_median
            df.cc_upper = intensity_result_full.cc_upper
            df.cc_se = intensity_result_full.cc_se
            df.ac1_lower = intensity_result_full.ac1_lower
            df.ac1_median = intensity_result_full.ac1_median
            df.ac1_upper = intensity_result_full.ac1_upper
            df.ac1_se = intensity_result_full.ac1_se
            df.ac2_lower = intensity_result_full.ac2_lower
            df.ac2_median = intensity_result_full.ac2_median
            df.ac2_upper = intensity_result_full.ac2_upper
            df.ac2_se = intensity_result_full.ac2_se
        end
    end

    # Add CorrelationTrait metadata columns (as scalars, repeated for each lag for CSV compatibility)
    # These specify how the empirical correlation functions were computed
    df.centering = fill(string(correlation_algorithm.centering), n_lags)
    df.multitau = fill(string(correlation_algorithm.multitau), n_lags)
    df.normalization = fill(string(correlation_algorithm.normalization), n_lags)
    df.biased = fill(correlation_algorithm.biased, n_lags)
    if correlation_algorithm.multitau == :multitau
        df.m = fill(correlation_algorithm.m, n_lags)
    end
    
    # Write to CSV (DataFrame constructor with kwargs preserves column order)
    CSV.write(output_file, df)
    println("Wrote empirical correlation functions to: $output_file")
end


# ============================================================================
# CORE CORRELATION FUNCTION COMPUTATION ENGINE
# ============================================================================
# Functions: compute_correlation_functions (main engine)
#            _validate_lags_for_traces, _compute_per_trace_correlation_functions, 
#            bootstrap_tracewise (helpers)
# Creation dates: compute_correlation_functions - refactored 2026-01-XX

# _validate_lags_for_traces moved to utilities.jl

"""
    compute_correlation_functions(traces, lags; mean1=nothing, mean2=nothing, correlation_algorithm=StandardCorrelation(), bootstrap=false, n_bootstrap=1000, frame_interval=1.0)

Compute cross-correlation and auto-correlation functions from a set of traces.

# Arguments
- `traces::Vector{Matrix{Float64}}`: Vector of trace matrices, each with 2 columns (x, y)
- `lags::Vector{<:Real}`: Vector of lags to compute (can be integers or floats)

# Keyword Arguments
- `mean1`: Global mean for first unit (if provided, used for centering when `correlation_algorithm.centering=:global_mean`)
- `mean2`: Global mean for second unit (if provided, used for centering when `correlation_algorithm.centering=:global_mean`)
- `correlation_algorithm::CorrelationTrait`: Algorithm traits (centering, multi-tau, normalization)
- `bootstrap::Bool`: Whether to compute bootstrap confidence intervals
- `n_bootstrap::Int`: Number of bootstrap iterations
- `frame_interval::Float64`: Sampling interval (default: 1.0)

# Returns
NamedTuple with:
- `cc`: Cross-correlation function C_XY(τ) for all lags
- `ac1`: Auto-correlation function C_XX(τ) for all lags
- `ac2`: Auto-correlation function C_YY(τ) for all lags
- `mean1`, `mean2`: Overall means (from provided means or computed from traces)
- `v1`, `v2`: Overall variances
- `lags`: The input lags vector
- Bootstrap fields (`*_lower`, `*_median`, `*_upper`, `*_se`) if `bootstrap=true`

# Notes
- If means are provided, they are used for centering when `correlation_algorithm.centering=:global_mean`
- Otherwise, means are computed from the traces
"""
function compute_correlation_functions(traces::Vector{Matrix{Float64}}, lags::Vector{<:Real}; mean1=nothing, mean2=nothing, correlation_algorithm=StandardCorrelation(), bootstrap::Bool=false, n_bootstrap::Int=1000, frame_interval::Float64=1.0)
    # Thin wrapper: dispatch to utilities.jl for all computation
    return compute_correlation_functions_traces(traces, lags; mean1=mean1, mean2=mean2, correlation_algorithm=correlation_algorithm, bootstrap=bootstrap, n_bootstrap=n_bootstrap, frame_interval=frame_interval)
end

# Helper functions moved to utilities.jl:
# - _validate_lags_for_traces
# - _prepare_centering_means  
# - _compute_windowed_means (replaced by compute_windowed_means)
# - _compute_per_trace_correlation_functions


# bootstrap_tracewise removed - use bootstrap_correlation_functions from utilities.jl instead
# This function has been replaced by bootstrap_correlation_function and bootstrap_correlation_functions in utilities.jl

# ############
# function compute_conditional_prob(traces::Vector{Matrix{Float64}}, lags::Vector{<:Real}; frame_interval=1.0)
#     n_traces = length(traces)
#     n_lags = length(lags)

#     # Numerator: Sum of X(t)Y(t+tau) across all traces
#     # Denominator: Sum of X(t) across all traces (weighted by available window T-tau)
#     numerator = zeros(n_lags)
#     denominator = zeros(n_lags)

#     for t_mat in traces
#         x = t_mat[:, 1]
#         y = t_mat[:, 2]
#         T = length(x)

#         for (li, lag) in enumerate(lags)
#             lag_frames = round(Int, lag / frame_interval)

#             # For a specific lag, find valid indices
#             if lag_frames >= 0
#                 idx_x = 1:(T - lag_frames)
#                 idx_y = (1 + lag_frames):T
#             else
#                 idx_x = (1 - lag_frames):T
#                 idx_y = 1:(T + lag_frames)
#             end

#             # Increment the conditional counts
#             # Numerator: total times both were ON
#             numerator[li] += sum(x[idx_x] .* y[idx_y])
#             # Denominator: total times the "source" (X) was ON
#             denominator[li] += sum(x[idx_x])
#         end
#     end

#     # The conditional probability P(Y=1 | X=1)
#     return numerator ./ denominator
# end
########################################################

# ============================================================================
# TRACE PREPARATION FUNCTIONS
# ============================================================================
# Functions: load_predicted_traces_csv (LEGACY), load_raw_reporters_csv, build_on_traces, 
#            prepare_traces_from_prediction
# Creation dates: 
#   - prepare_traces_from_prediction, load_raw_reporters_csv, build_on_traces: 2026-01-06 (AI-created)
#   - load_predicted_traces_csv: OLD (pre-2025, legacy, kept for compatibility - may be deprecated)
# Purpose: Load and prepare traces from CSV files for covariance computation
# NOTE: prepare_traces_from_prediction is the canonical entrypoint (recommended)
# NOTE: load_predicted_traces_csv is legacy and may be deprecated

"""
    load_predicted_traces_csv(unit1_filename, unit2_filename; reporters_col=nothing, trace_id_col=nothing)

Find a rate file for a specific gene from a list of files.

# Arguments
- `files::Vector{String}`: List of file names to search
- `gene::String`: Gene name to search for

# Returns
- `String` or `Int`: File name if found, 0 if not found

# Notes
- Searches for files containing "_gene_" pattern
- Returns the first matching file
- Returns 0 if no matching file is found
- Assumes gene names are embedded in filenames with underscores
- Useful for finding rate files for specific genes

# Examples
```julia
# Find rate file for MYC gene
files = ["rates_MYC_control.txt", "rates_FOS_treatment.txt"]
rate_file = getratefile(files, "MYC")  # Returns: "rates_MYC_control.txt"

# Gene not found
rate_file = getratefile(files, "UNKNOWN")  # Returns: 0
```
"""


"""
    getratefile(folder, G, cond)

Find rate files in a folder for a specific model and condition.

# Arguments
- `folder::String`: Path to folder to search
- `G`: Model identifier
- `cond::String`: Condition identifier

# Returns
- `Vector{String}`: Vector of matching file names

# Notes
- Searches for files containing "rates_", condition, and model identifiers
- Filters files by multiple criteria: rates prefix, condition, and model
- Returns all files that match all criteria
- Assumes files follow standard naming conventions
- Useful for finding all rate files for a specific model-condition combination

# Examples
```julia
# Find rate files for G=2 model in control condition
files = getratefile("results/", "2", "control")

# Returns vector of files like: ["rates_gene1_control_2.txt", "rates_gene2_control_2.txt"]
```
"""

"""
    change_name(folder, oldname, newname)

Rename all files in a folder that contain a specific pattern.

# Arguments
- `folder::String`: Path to folder containing files
- `oldname::String`: Old name pattern to replace
- `newname::String`: New name pattern to use

# Returns
- `Nothing`: Renames files in place

# Notes
- Searches for files containing the oldname pattern
- Replaces oldname with newname in all matching filenames
- Uses force=true to overwrite existing files if necessary
- Renames files in the same folder
- Useful for batch renaming files with consistent patterns

# Examples
```julia
# Rename all files containing "old_condition" to "new_condition"
change_name("results/", "old_condition", "new_condition")

# This would rename:
# - "rates_gene1_old_condition_2.txt" → "rates_gene1_new_condition_2.txt"
# - "stats_gene2_old_condition_2.txt" → "stats_gene2_new_condition_2.txt"
```
"""

"""
    get_histogram_rna(gene, datacond, datapath)

Get normalized RNA histogram for a specific gene and condition.

# Arguments
- `gene::String`: Gene name
- `datacond::String`: Data condition identifier
- `datapath::String`: Path to data directory

# Returns
- `Vector{Float64}`: Normalized RNA histogram

# Notes
- Reads RNA data for the specified gene and condition
- Normalizes the histogram to sum to 1
- Uses read_rna function to load data
- Uses normalize_histogram function for normalization
- Returns probability distribution over RNA counts

# Examples
```julia
# Get normalized RNA histogram for MYC gene in control condition
hist = get_histogram_rna("MYC", "control", "data/")

# hist contains normalized probabilities for each RNA count
```
"""

"""
    make_vector(x, n)

Create a vector of length n, either by repeating a scalar or using an existing vector.

# Arguments
- `x`: Input value (scalar or vector)
- `n::Int`: Desired length of output vector

# Returns
- `Vector`: Vector of length n

# Notes
- If x is a scalar: creates a vector of length n filled with x
- If x is already a vector: returns x as-is
- Useful for ensuring consistent vector lengths in function calls
- Handles both numeric and non-numeric types

# Examples
```julia
# Create vector from scalar
make_vector(5, 3)  # Returns: [5, 5, 5]

# Pass through existing vector
make_vector([1, 2, 3], 3)  # Returns: [1, 2, 3]

# Create vector of strings
make_vector("default", 2)  # Returns: ["default", "default"]
```
"""

"""
    separate_dataframe(df, G)

Separate a DataFrame into multiple DataFrames based on gene states.

# Arguments
- `df::DataFrame`: Input DataFrame
- `G::Int`: Number of gene states

# Returns
- `Vector{DataFrame}`: Vector of separated DataFrames

# Notes
- Creates separate DataFrame for each gene state
- Extracts appropriate columns for each state based on G
- Assumes column structure follows standard format
- Uses nsets variable (should be defined in scope)
- Simplified version without condition handling

# Examples
```julia
# Separate DataFrame by gene states
dfs = separate_dataframe(df, 3)

# Returns vector of DataFrames, one for each gene state
```
"""

"""
    compute_deviance(outfile, ratefile::String, cond, n, datapath, root)

Compute deviance for all genes in a rate file and write results to output file.

# Arguments
- `outfile::String`: Path to output file
- `ratefile::String`: Path to rate file
- `cond::String`: Condition identifier
- `n::Int`: Number of states
- `datapath::String`: Path to data directory
- `root::String`: Root path for data

# Returns
- `Nothing`: Writes deviance results to outfile

# Notes
- Reads rate file and computes deviance for each gene
- Uses deviance function for individual gene calculations
- Writes results in format: gene_name, deviance_value
- Assumes rate file has gene names in first column
- Useful for batch deviance computation across multiple genes

# Examples
```julia
# Compute deviance for all genes in rate file
compute_deviance("deviance_results.txt", "rates.txt", "control", 2, "data/", ".")

# Creates output file with deviance values for each gene
```
"""

"""
    data_covariance(traces, lags)

Compute autocovariance and cross-covariance from trace data.

# Arguments
- `traces::Vector`: Vector of trace data, where each trace is a matrix with columns [enhancer, gene]
- `lags::Vector{Int}`: Vector of time lags

# Returns
- `Tuple{Vector, Vector, Vector, Vector}`: (ac1, ac2, cov, lags)
  - `ac1`: Autocovariance function for first trace (enhancer)
  - `ac2`: Autocovariance function for second trace (gene)
  - `cov`: Cross-covariance function between traces
  - `lags`: Time lags (same as input)

# Notes
- Computes autocovariance for first and second traces
- Computes cross-covariance between traces
- Averages results across all traces
- Uses `crosscorrelation_function` with empirical means which:
  - Removes mean (handles trends/drift)
  - Only uses valid pairs (handles edge effects automatically)
  - For lag τ, only uses pairs (t, t+τ) where both indices are valid
- Edge effects are automatically handled by only computing covariances
  for valid time pairs within each trace

# Examples
```julia
# Compute covariance functions from trace data
lags = collect(-60:1:60)
ac1, ac2, cov, lags = data_covariance(traces, lags)

# ac1, ac2: autocovariance functions
# cov: cross-covariance function
```
"""


"""
    load_predicted_traces_csv(filename; enhancer_col="Reporters1_n", gene_col="Reporters2_n", trace_id_col=nothing)

Load predicted traces from CSV file and convert reporter counts to binary ON/OFF states.

# Arguments
- `filename::String`: Path to CSV file
- `enhancer_col::String="Reporters1_n"`: Column name for enhancer reporter counts
- `gene_col::String="Reporters2_n"`: Column name for gene reporter counts
- `trace_id_col::Union{String, Nothing}=nothing`: Column name for trace ID (if multiple traces in file)

# Returns
- `Vector{Matrix{Float64}}`: Vector of binary ON/OFF trace matrices, each with columns [enhancer_ONOFF, gene_ONOFF]
  - ON = 1 (reporters > 0), OFF = 0 (reporters <= 0)

# Notes
- If `trace_id_col` is provided, traces are grouped by that column
- If not provided, assumes single trace or uses row index to group traces
- Binary conversion: > 0 reporter = ON (1), <= 0 reporter = OFF (0)
"""
function load_predicted_traces_csv(unit1_filename, unit2_filename; reporters_col=nothing, trace_id_col=nothing)
    # Load CSV files (requires CSV and DataFrames packages)
    # unit1_filename: file for enhancer (unit 1)
    # unit2_filename: file for gene (unit 2)
    df1 = CSV.read(unit1_filename, DataFrame)
    df2 = CSV.read(unit2_filename, DataFrame)

    # Auto-detect reporter columns if not specified
    if isnothing(reporters_col)
        # Look for columns matching pattern "Reporters_1", "Reporters_2", etc.
        reporters_cols1 = filter(n -> startswith(string(n), "Reporters_"), names(df1))
        reporters_cols2 = filter(n -> startswith(string(n), "Reporters_"), names(df2))

        if isempty(reporters_cols1) || isempty(reporters_cols2)
            error("Could not find 'Reporters_n' columns in files. Found columns: unit1=$(names(df1)), unit2=$(names(df2))")
        end

        # Sort to ensure matching order (Reporters_1, Reporters_2, ...)
        reporters_cols1 = sort(reporters_cols1)
        reporters_cols2 = sort(reporters_cols2)

        if length(reporters_cols1) != length(reporters_cols2)
            error("Different number of reporter columns: unit1 has $(length(reporters_cols1)), unit2 has $(length(reporters_cols2))")
        end

        # Extract all traces
        traces = Matrix{Float64}[]
        for i in 1:length(reporters_cols1)
            unit1_counts = df1[:, reporters_cols1[i]]
            unit2_counts = df2[:, reporters_cols2[i]]

            # Handle missing values: replace with 0 (OFF state)
            unit1_counts = coalesce.(unit1_counts, 0.0)
            unit2_counts = coalesce.(unit2_counts, 0.0)

            # Convert to binary ON/OFF: > 0 = ON (1), <= 0 = OFF (0)
            # Match experimentalist convention: reporter > 0
            unit1_onoff = Float64.(unit1_counts .> 0)
            unit2_onoff = Float64.(unit2_counts .> 0)

            # Ensure same length
            if length(unit1_onoff) != length(unit2_onoff)
                error("Trace $i has different lengths: unit1=$(length(unit1_onoff)), unit2=$(length(unit2_onoff))")
            end

            trace_matrix = hcat(unit1_onoff, unit2_onoff)
            push!(traces, trace_matrix)
        end

        return traces
    end

    # Original behavior: single reporter column specified
    unit1_counts = df1[:, reporters_col]
    unit2_counts = df2[:, reporters_col]

    # Handle missing values: replace with 0 (OFF state)
    unit1_counts = coalesce.(unit1_counts, 0.0)
    unit2_counts = coalesce.(unit2_counts, 0.0)

    # Convert to binary ON/OFF: > 0 = ON (1), <= 0 = OFF (0)
    # Match experimentalist convention: reporter > 0
    unit1_onoff = Float64.(unit1_counts .> 0)
    unit2_onoff = Float64.(unit2_counts .> 0)

    # Group by trace if trace_id_col is provided
    if !isnothing(trace_id_col)
        # Check if trace_id_col exists in both dataframes
        if trace_id_col in names(df1) && trace_id_col in names(df2)
            trace_ids1 = df1[:, trace_id_col]
            trace_ids2 = df2[:, trace_id_col]
            unique_traces1 = unique(trace_ids1)
            unique_traces2 = unique(trace_ids2)

            # Ensure both files have the same trace IDs
            if unique_traces1 != unique_traces2
                error("Trace IDs in unit1 and unit2 files do not match")
            end

            traces = Matrix{Float64}[]
            for trace_id in unique_traces1
                idx1 = trace_ids1 .== trace_id
                idx2 = trace_ids2 .== trace_id
                # Ensure same length for this trace
                if sum(idx1) != sum(idx2)
                    error("Trace $trace_id has different lengths in unit1 and unit2 files")
                end
                trace_matrix = hcat(unit1_onoff[idx1], unit2_onoff[idx2])
                push!(traces, trace_matrix)
            end

            return traces
        else
            error("trace_id_col '$trace_id_col' not found in one or both dataframes")
        end
    else
        # Assume single trace - ensure both files have same number of rows
        if length(unit1_onoff) != length(unit2_onoff)
            error("Unit1 and unit2 files have different numbers of rows")
        end
        trace_matrix = hcat(unit1_onoff, unit2_onoff)
        return [trace_matrix]
    end
end

"""
    load_raw_reporters_csv(unit1_filename, unit2_filename; reporters_col=nothing, trace_id_col=nothing)

Load raw reporter counts from CSV files (without converting to binary ON/OFF).

# Arguments
- `unit1_filename::String`: Path to CSV file for unit 1 (enhancer)
- `unit2_filename::String`: Path to CSV file for unit 2 (gene)
- `reporters_col::Union{String, Nothing}=nothing`: Column name for reporter counts (auto-detected if nothing)
- `trace_id_col::Union{String, Nothing}=nothing`: Column name for trace ID (if multiple traces)

# Returns
- `Vector{Matrix{Float64}}`: Vector of raw reporter count matrices, each with columns [unit1_counts, unit2_counts]
"""
function load_raw_reporters_csv(unit1_filename, unit2_filename; reporters_col=nothing, trace_id_col=nothing)
    df1 = CSV.read(unit1_filename, DataFrame)
    df2 = CSV.read(unit2_filename, DataFrame)

    # Auto-detect reporter columns if not specified
    if isnothing(reporters_col)
        reporters_cols1 = filter(n -> startswith(string(n), "Reporters_"), names(df1))
        reporters_cols2 = filter(n -> startswith(string(n), "Reporters_"), names(df2))

        if isempty(reporters_cols1) || isempty(reporters_cols2)
            error("Could not find 'Reporters_' columns in files. Found columns: unit1=$(names(df1)), unit2=$(names(df2))")
        end

        reporters_cols1 = sort(reporters_cols1)
        reporters_cols2 = sort(reporters_cols2)

        if length(reporters_cols1) != length(reporters_cols2)
            error("Different number of reporter columns: unit1 has $(length(reporters_cols1)), unit2 has $(length(reporters_cols2))")
        end

        traces = Matrix{Float64}[]
        for i in 1:length(reporters_cols1)
            unit1_counts = df1[:, reporters_cols1[i]]
            unit2_counts = df2[:, reporters_cols2[i]]

            # Handle missing values: replace with 0
            unit1_counts = coalesce.(unit1_counts, 0.0)
            unit2_counts = coalesce.(unit2_counts, 0.0)

            # Keep as raw counts (don't convert to binary)
            if length(unit1_counts) != length(unit2_counts)
                error("Trace $i has different lengths: unit1=$(length(unit1_counts)), unit2=$(length(unit2_counts))")
            end

            trace_matrix = hcat(unit1_counts, unit2_counts)
            push!(traces, trace_matrix)
        end

        return traces
    end

    # Single reporter column specified
    unit1_counts = df1[:, reporters_col]
    unit2_counts = df2[:, reporters_col]

    # Handle missing values
    unit1_counts = coalesce.(unit1_counts, 0.0)
    unit2_counts = coalesce.(unit2_counts, 0.0)

    # Group by trace if trace_id_col is provided
    if !isnothing(trace_id_col)
        if trace_id_col in names(df1) && trace_id_col in names(df2)
            trace_ids1 = df1[:, trace_id_col]
            trace_ids2 = df2[:, trace_id_col]
            unique_traces1 = unique(trace_ids1)
            unique_traces2 = unique(trace_ids2)

            if unique_traces1 != unique_traces2
                error("Trace IDs in unit1 and unit2 files do not match")
            end

            traces = Matrix{Float64}[]
            for trace_id in unique_traces1
                idx1 = trace_ids1 .== trace_id
                idx2 = trace_ids2 .== trace_id
                if sum(idx1) != sum(idx2)
                    error("Trace $trace_id has different lengths in unit1 and unit2 files")
                end
                trace_matrix = hcat(unit1_counts[idx1], unit2_counts[idx2])
                push!(traces, trace_matrix)
            end

            return traces
        else
            error("trace_id_col '$trace_id_col' not found in one or both dataframes")
        end
    else
        # Single trace
        if length(unit1_counts) != length(unit2_counts)
            error("Unit1 and unit2 files have different numbers of rows")
        end
        trace_matrix = hcat(unit1_counts, unit2_counts)
        return [trace_matrix]
    end
end

"""
    build_on_traces(reporter_traces; threshold=0)

Convert reporter-count traces to binary ON/OFF traces.

Each trace is a matrix with columns `[unit1, unit2]` (reporter counts). The returned traces
are Float64 matrices with the same shape, where ON is defined as `count > threshold`.
"""
function build_on_traces(reporter_traces::Vector{<:AbstractMatrix}; threshold=0)
    on_traces = Matrix{Float64}[]
    for t in reporter_traces
        # treat missing as 0 (OFF)
        r1 = coalesce.(t[:, 1], 0.0)
        r2 = coalesce.(t[:, 2], 0.0)
        push!(on_traces, hcat(Float64.(r1 .> threshold), Float64.(r2 .> threshold)))
    end
    return on_traces
end

"""
    prepare_traces_from_prediction(unit1_filename, unit2_filename; reporters_col=nothing, trace_id_col=nothing, on_threshold=0)

Prepare a *uniform* trace bundle from prediction CSV files.

Returns both:
- `reporter_traces`: raw reporter counts (Float64)
- `on_traces`: binary ON/OFF derived from reporter counts (`count > on_threshold`)

This is the canonical entrypoint for trace preparation from the prediction CSV format.
"""
function prepare_traces_from_prediction(unit1_filename, unit2_filename; reporters_col=nothing, trace_id_col=nothing, on_threshold=0)
    reporter_traces = load_raw_reporters_csv(unit1_filename, unit2_filename; reporters_col=reporters_col, trace_id_col=trace_id_col)
    on_traces = build_on_traces(reporter_traces; threshold=on_threshold)
    return (reporter_traces=reporter_traces, on_traces=on_traces)
end


# ============================================================================
# MODEL SCORING FUNCTIONS
# ============================================================================
# Functions: score_models_from_traces, compute_score_metrics, summarize_model_scores
# Creation dates:
#   - score_models_from_traces: 2025-12-05 or earlier (user-created, exported)
#   - compute_score_metrics: 2026-01-06 (AI-created helper)
#   - summarize_model_scores: OLD (pre-2025, user-created)
# Purpose: Score empirical vs theoretical covariance functions
# NOTE: score_models_from_traces is the main user-facing function (exported)

"""
    score_models_from_traces(enhancer_file, gene_file, crosscov_folder; 
        crosscov_pattern="crosscovariance_tracejoint-HBEC-nstate_enhancer-gene")

Score model predictions against empirical cross-covariance by reading precomputed results from CSV files.

The crosscovariance files contain the theoretical predictions (precomputed by `write_xcorr`).
The empirical file contains empirical results (precomputed by `write_xcorr_empirical`).
Workflow:
1. Read empirical cross-covariance from CSV file (precomputed)
2. Read theoretical cross-covariance from CSV files in crosscov_folder (includes tau/lags)
3. Score theoretical predictions (from files) against empirical data (from file)

The coupling model is extracted from crosscovariance filenames (two characters after "gene").
Lags are read from the `tau` column in the crosscovariance files.

# Arguments
- `empirical_file::String`: Path to empirical cross-covariance CSV file (from `write_xcorr_empirical`)
- `crosscov_folder::String`: Folder containing crosscovariance CSV files (with theoretical predictions from `write_xcorr`)
- `crosscov_pattern::String`: Filename pattern for crosscovariance files

# Returns
Dictionary mapping coupling model identifiers to named tuples with:
- `coupling_model::String`: Coupling model identifier
- `cc_empirical`: Empirical cross-covariance computed from traces
- `cc_theory`: Theoretical cross-covariance read from file
- `cc_l2_norm`: L² norm (theory vs empirical)
- `cc_linf_norm`: L∞ norm (theory vs empirical)
- `mean1_empirical`, `mean2_empirical`: Empirical ON state means
- `lags`: Time lags

# Example
```julia
results = score_models_from_traces(
    "/path/to/enhancer/folder",
    "/path/to/gene/folder",
    r, transitions, G, R, S, insertstep, interval, probfn, coupling, lags
)

# Access results for coupling model "31":
results["31"].cc_l2_norm
results["31"].cc_empirical
results["31"].cc_theory
```
"""
function score_models_from_traces(empirical_file::String, crosscov_folder::String;
    crosscov_pattern="crosscorrelation_tracejoint-HBEC-nstate_enhancer-gene")

    # Find all crosscorrelation files (these contain the coupling model identifier)
    crosscov_files = filter(f -> startswith(f, crosscov_pattern) && endswith(f, ".csv"), readdir(crosscov_folder))

    if isempty(crosscov_files)
        error("No crosscorrelation files matching pattern '$crosscov_pattern*.csv' found in folder '$crosscov_folder'")
    end

    # Read empirical results from CSV file (precomputed by write_correlation_functions_empirical)
    df_empirical = CSV.read(empirical_file, DataFrame)
    if !("tau" in names(df_empirical))
        error("Empirical file missing 'tau' column")
    end
    lags_empirical_full = Vector{Float64}(df_empirical.tau)
    
    # Read CorrelationTrait metadata from empirical file (if present)
    # Defaults to StandardCorrelation() if not found (backward compatibility)
    if hasproperty(df_empirical, :centering) && hasproperty(df_empirical, :multitau) && hasproperty(df_empirical, :normalization) && hasproperty(df_empirical, :biased)
        empirical_centering = Symbol(df_empirical.centering[1])
        empirical_multitau = Symbol(df_empirical.multitau[1])
        empirical_normalization = Symbol(df_empirical.normalization[1])
        empirical_biased = df_empirical.biased[1]
        empirical_m = hasproperty(df_empirical, :m) ? df_empirical.m[1] : 16
        empirical_correlation_algorithm = CorrelationTrait(
            centering=empirical_centering,
            multitau=empirical_multitau,
            normalization=empirical_normalization,
            biased=empirical_biased,
            m=empirical_m
        )
    else
        # Backward compatibility: use StandardCorrelation() if metadata not present
        @warn "Empirical file missing CorrelationTrait metadata. Using StandardCorrelation() as default. Consider regenerating empirical file with updated code."
        empirical_correlation_algorithm = StandardCorrelation()
    end

    # Extract empirical ON state data (will be filtered per file to match theory range)
    data_result_full = (
        cc_unnormalized=df_empirical.cc_ON,
        ac1_unnormalized=df_empirical.ac1_ON,
        ac2_unnormalized=df_empirical.ac2_ON,
        mean1=df_empirical.m_ON1[1],
        mean2=df_empirical.m_ON2[1],
        lags=lags_empirical_full,
        cc_unnormalized_lower=hasproperty(df_empirical, :cc_ON_lower) ? df_empirical.cc_ON_lower : nothing,
        cc_unnormalized_median=hasproperty(df_empirical, :cc_ON_median) ? df_empirical.cc_ON_median : nothing,
        cc_unnormalized_upper=hasproperty(df_empirical, :cc_ON_upper) ? df_empirical.cc_ON_upper : nothing,
        cc_unnormalized_se=hasproperty(df_empirical, :cc_ON_se) ? df_empirical.cc_ON_se : nothing,
        ac1_unnormalized_lower=hasproperty(df_empirical, :ac1_ON_lower) ? df_empirical.ac1_ON_lower : nothing,
        ac1_unnormalized_median=hasproperty(df_empirical, :ac1_ON_median) ? df_empirical.ac1_ON_median : nothing,
        ac1_unnormalized_upper=hasproperty(df_empirical, :ac1_ON_upper) ? df_empirical.ac1_ON_upper : nothing,
        ac1_unnormalized_se=hasproperty(df_empirical, :ac1_ON_se) ? df_empirical.ac1_ON_se : nothing,
        ac2_unnormalized_lower=hasproperty(df_empirical, :ac2_ON_lower) ? df_empirical.ac2_ON_lower : nothing,
        ac2_unnormalized_median=hasproperty(df_empirical, :ac2_ON_median) ? df_empirical.ac2_ON_median : nothing,
        ac2_unnormalized_upper=hasproperty(df_empirical, :ac2_ON_upper) ? df_empirical.ac2_ON_upper : nothing,
        ac2_unnormalized_se=hasproperty(df_empirical, :ac2_ON_se) ? df_empirical.ac2_ON_se : nothing
    )

    # Extract empirical Reporter data (will be filtered per file to match theory range)
    reporter_result_full = (
        cc=hasproperty(df_empirical, :cc_Reporters) ? df_empirical.cc_Reporters : nothing,
        ac1=hasproperty(df_empirical, :ac1_Reporters) ? df_empirical.ac1_Reporters : nothing,
        ac2=hasproperty(df_empirical, :ac2_Reporters) ? df_empirical.ac2_Reporters : nothing,
        mean1=hasproperty(df_empirical, :m_Reporters1) ? df_empirical.m_Reporters1[1] : nothing,
        mean2=hasproperty(df_empirical, :mReporters2) ? df_empirical.mReporters2[1] : nothing,
        lags=lags_empirical_full,
        cc_lower=hasproperty(df_empirical, :cc_Reporters_lower) ? df_empirical.cc_Reporters_lower : nothing,
        cc_median=hasproperty(df_empirical, :cc_Reporters_median) ? df_empirical.cc_Reporters_median : nothing,
        cc_upper=hasproperty(df_empirical, :cc_Reporters_upper) ? df_empirical.cc_Reporters_upper : nothing,
        cc_se=hasproperty(df_empirical, :cc_Reporters_se) ? df_empirical.cc_Reporters_se : nothing,
        ac1_lower=hasproperty(df_empirical, :ac1_Reporters_lower) ? df_empirical.ac1_Reporters_lower : nothing,
        ac1_median=hasproperty(df_empirical, :ac1_Reporters_median) ? df_empirical.ac1_Reporters_median : nothing,
        ac1_upper=hasproperty(df_empirical, :ac1_Reporters_upper) ? df_empirical.ac1_Reporters_upper : nothing,
        ac1_se=hasproperty(df_empirical, :ac1_Reporters_se) ? df_empirical.ac1_Reporters_se : nothing,
        ac2_lower=hasproperty(df_empirical, :ac2_Reporters_lower) ? df_empirical.ac2_Reporters_lower : nothing,
        ac2_median=hasproperty(df_empirical, :ac2_Reporters_median) ? df_empirical.ac2_Reporters_median : nothing,
        ac2_upper=hasproperty(df_empirical, :ac2_Reporters_upper) ? df_empirical.ac2_Reporters_upper : nothing,
        ac2_se=hasproperty(df_empirical, :ac2_Reporters_se) ? df_empirical.ac2_Reporters_se : nothing
    )

    results = Dict{String,NamedTuple}()

    for crosscov_file in crosscov_files
        # Extract coupling model from crosscovariance filename
        gene_pos = findfirst("gene", crosscov_file)
        if isnothing(gene_pos)
            @warn "Could not find 'gene' in filename '$crosscov_file', skipping"
            continue
        end

        coupling_start = gene_pos.stop + 1
        if coupling_start + 1 > length(crosscov_file)
            @warn "Filename '$crosscov_file' too short to extract coupling model, skipping"
            continue
        end

        coupling_model_str = crosscov_file[coupling_start:coupling_start+1]
        # Parse model string (G,R,S,insertstep) from filename: ..._gene_MYC_model_nalleles.csv
        base_no_ext = endswith(crosscov_file, ".csv") ? crosscov_file[1:end-4] : crosscov_file
        v = split(base_no_ext, "_")
        model_str = length(v) >= 6 ? v[5] : ""
        results_key = isempty(model_str) ? coupling_model_str : coupling_model_str * "_" * model_str
        crosscov_filepath = joinpath(crosscov_folder, crosscov_file)

        # 2. Read crosscovariance from file
        try
            df = CSV.read(crosscov_filepath, DataFrame)
            df_names_str = string.(names(df))

            if !("tau" in df_names_str)
                @warn "File '$crosscov_file' missing required column 'tau', skipping"
                continue
            end

            tau_data = Vector{Float64}(df.tau)
            
            # Use shortest overlapping lag range (intersection of empirical and theory ranges)
            # This ensures we compare theory and empirical on the same lags
            tau_min_empirical = minimum(lags_empirical_full)
            tau_max_empirical = maximum(lags_empirical_full)
            tau_min_theory = minimum(tau_data)
            tau_max_theory = maximum(tau_data)
            
            tau_min = max(tau_min_empirical, tau_min_theory)
            tau_max = min(tau_max_empirical, tau_max_theory)
            
            if tau_min > tau_max
                @warn "No overlapping lags between empirical (range: $tau_min_empirical to $tau_max_empirical) and theory (range: $tau_min_theory to $tau_max_theory) in file '$crosscov_file', skipping"
                continue
            end
            
            # Filter both empirical and theory lags to the overlapping range
            lags_empirical_filtered = filter(lag -> tau_min <= lag <= tau_max, lags_empirical_full)
            lags_theory_filtered = filter(lag -> tau_min <= lag <= tau_max, tau_data)
            
            # Use the empirical lags as the reference (they define the sampling)
            lags = lags_empirical_filtered
            
            # Filter empirical data to match filtered lags
            empirical_mask = [tau_min <= lag <= tau_max for lag in lags_empirical_full]
            
            # Create filtered empirical data structures for this theory file
            data_result = (
                cc_unnormalized=Vector{Float64}(data_result_full.cc_unnormalized[empirical_mask]),
                ac1_unnormalized=Vector{Float64}(data_result_full.ac1_unnormalized[empirical_mask]),
                ac2_unnormalized=Vector{Float64}(data_result_full.ac2_unnormalized[empirical_mask]),
                mean1=data_result_full.mean1,
                mean2=data_result_full.mean2,
                lags=lags,
                cc_unnormalized_lower=hasproperty(data_result_full, :cc_unnormalized_lower) && !isnothing(data_result_full.cc_unnormalized_lower) ? Vector{Float64}(data_result_full.cc_unnormalized_lower[empirical_mask]) : nothing,
                cc_unnormalized_median=hasproperty(data_result_full, :cc_unnormalized_median) && !isnothing(data_result_full.cc_unnormalized_median) ? Vector{Float64}(data_result_full.cc_unnormalized_median[empirical_mask]) : nothing,
                cc_unnormalized_upper=hasproperty(data_result_full, :cc_unnormalized_upper) && !isnothing(data_result_full.cc_unnormalized_upper) ? Vector{Float64}(data_result_full.cc_unnormalized_upper[empirical_mask]) : nothing,
                cc_unnormalized_se=hasproperty(data_result_full, :cc_unnormalized_se) && !isnothing(data_result_full.cc_unnormalized_se) ? Vector{Float64}(data_result_full.cc_unnormalized_se[empirical_mask]) : nothing,
                ac1_unnormalized_lower=hasproperty(data_result_full, :ac1_unnormalized_lower) && !isnothing(data_result_full.ac1_unnormalized_lower) ? Vector{Float64}(data_result_full.ac1_unnormalized_lower[empirical_mask]) : nothing,
                ac1_unnormalized_median=hasproperty(data_result_full, :ac1_unnormalized_median) && !isnothing(data_result_full.ac1_unnormalized_median) ? Vector{Float64}(data_result_full.ac1_unnormalized_median[empirical_mask]) : nothing,
                ac1_unnormalized_upper=hasproperty(data_result_full, :ac1_unnormalized_upper) && !isnothing(data_result_full.ac1_unnormalized_upper) ? Vector{Float64}(data_result_full.ac1_unnormalized_upper[empirical_mask]) : nothing,
                ac1_unnormalized_se=hasproperty(data_result_full, :ac1_unnormalized_se) && !isnothing(data_result_full.ac1_unnormalized_se) ? Vector{Float64}(data_result_full.ac1_unnormalized_se[empirical_mask]) : nothing,
                ac2_unnormalized_lower=hasproperty(data_result_full, :ac2_unnormalized_lower) && !isnothing(data_result_full.ac2_unnormalized_lower) ? Vector{Float64}(data_result_full.ac2_unnormalized_lower[empirical_mask]) : nothing,
                ac2_unnormalized_median=hasproperty(data_result_full, :ac2_unnormalized_median) && !isnothing(data_result_full.ac2_unnormalized_median) ? Vector{Float64}(data_result_full.ac2_unnormalized_median[empirical_mask]) : nothing,
                ac2_unnormalized_upper=hasproperty(data_result_full, :ac2_unnormalized_upper) && !isnothing(data_result_full.ac2_unnormalized_upper) ? Vector{Float64}(data_result_full.ac2_unnormalized_upper[empirical_mask]) : nothing,
                ac2_unnormalized_se=hasproperty(data_result_full, :ac2_unnormalized_se) && !isnothing(data_result_full.ac2_unnormalized_se) ? Vector{Float64}(data_result_full.ac2_unnormalized_se[empirical_mask]) : nothing
            )
            
            reporter_result = (
                cc=!isnothing(reporter_result_full.cc) ? Vector{Float64}(reporter_result_full.cc[empirical_mask]) : nothing,
                ac1=!isnothing(reporter_result_full.ac1) ? Vector{Float64}(reporter_result_full.ac1[empirical_mask]) : nothing,
                ac2=!isnothing(reporter_result_full.ac2) ? Vector{Float64}(reporter_result_full.ac2[empirical_mask]) : nothing,
                mean1=reporter_result_full.mean1,
                mean2=reporter_result_full.mean2,
                lags=lags,
                cc_lower=!isnothing(reporter_result_full.cc_lower) ? Vector{Float64}(reporter_result_full.cc_lower[empirical_mask]) : nothing,
                cc_median=!isnothing(reporter_result_full.cc_median) ? Vector{Float64}(reporter_result_full.cc_median[empirical_mask]) : nothing,
                cc_upper=!isnothing(reporter_result_full.cc_upper) ? Vector{Float64}(reporter_result_full.cc_upper[empirical_mask]) : nothing,
                cc_se=!isnothing(reporter_result_full.cc_se) ? Vector{Float64}(reporter_result_full.cc_se[empirical_mask]) : nothing,
                ac1_lower=!isnothing(reporter_result_full.ac1_lower) ? Vector{Float64}(reporter_result_full.ac1_lower[empirical_mask]) : nothing,
                ac1_median=!isnothing(reporter_result_full.ac1_median) ? Vector{Float64}(reporter_result_full.ac1_median[empirical_mask]) : nothing,
                ac1_upper=!isnothing(reporter_result_full.ac1_upper) ? Vector{Float64}(reporter_result_full.ac1_upper[empirical_mask]) : nothing,
                ac1_se=!isnothing(reporter_result_full.ac1_se) ? Vector{Float64}(reporter_result_full.ac1_se[empirical_mask]) : nothing,
                ac2_lower=!isnothing(reporter_result_full.ac2_lower) ? Vector{Float64}(reporter_result_full.ac2_lower[empirical_mask]) : nothing,
                ac2_median=!isnothing(reporter_result_full.ac2_median) ? Vector{Float64}(reporter_result_full.ac2_median[empirical_mask]) : nothing,
                ac2_upper=!isnothing(reporter_result_full.ac2_upper) ? Vector{Float64}(reporter_result_full.ac2_upper[empirical_mask]) : nothing,
                ac2_se=!isnothing(reporter_result_full.ac2_se) ? Vector{Float64}(reporter_result_full.ac2_se[empirical_mask]) : nothing
            )
            
            # Initialize all variables at the start to avoid scoping issues
            cc_ON_theory_unnorm_all = nothing
            mON1_theory = nothing
            mON2_theory = nothing
            ac1_ON_theory_all = nothing
            ac2_ON_theory_all = nothing
            cc_Reporters_theory_unnorm_all = nothing
            mR1_theory = nothing
            mR2_theory = nothing
            ac1_Reporters_theory_all = nothing
            ac2_Reporters_theory_all = nothing

            # Check for ON state columns
            has_cc_ON_unnorm = "cc_ON_unnormalized" in df_names_str
            has_m_ON1 = "m_ON1" in df_names_str
            has_m_ON2 = "m_ON2" in df_names_str
            has_cc_ON = "cc_ON" in df_names_str
            has_ac1_ON = "ac1_ON" in df_names_str
            has_ac2_ON = "ac2_ON" in df_names_str

            # Check for Reporter columns
            has_cc_Reporters_unnorm = "cc_Reporters_unnormalized" in df_names_str
            has_m_Reporters1 = "m_Reporters1" in df_names_str
            has_mReporters2 = "mReporters2" in df_names_str
            has_cc_Reporters = "cc_Reporters" in df_names_str
            has_ac1_Reporters = "ac1_Reporters" in df_names_str
            has_ac2_Reporters = "ac2_Reporters" in df_names_str

            if !has_cc_ON && !has_cc_ON_unnorm
                @warn "File '$crosscov_file' missing required columns for ON state ('cc_ON' or 'cc_ON_unnormalized'), skipping"
                continue
            end

            # Read ON state theoretical predictions (unnormalized)
            # Theory files contain raw uncentered unnormalized correlation functions (E[xy])
            if has_cc_ON_unnorm
                cc_ON_theory_unnorm_all = df.cc_ON_unnormalized
                mON1_theory = has_m_ON1 ? df.m_ON1[1] : nothing
                mON2_theory = has_m_ON2 ? df.m_ON2[1] : nothing
            elseif has_cc_ON
                # cc_ON from write_correlation_functions is raw uncentered (E[xy]), use it directly
                cc_ON_theory_unnorm_all = df.cc_ON
                mON1_theory = has_m_ON1 ? df.m_ON1[1] : nothing
                mON2_theory = has_m_ON2 ? df.m_ON2[1] : nothing
            else
                @warn "File '$crosscov_file' missing required ON state columns ('cc_ON' or 'cc_ON_unnormalized'), skipping"
                continue
            end

            # Read ON state autocovariances if available
            if has_ac1_ON
                ac1_ON_theory_all = df.ac1_ON
            end
            if has_ac2_ON
                ac2_ON_theory_all = df.ac2_ON
            end

            # Read Reporter theoretical predictions (unnormalized)
            # Theory files contain raw uncentered unnormalized correlation functions (E[xy])
            if has_cc_Reporters_unnorm
                cc_Reporters_theory_unnorm_all = df.cc_Reporters_unnormalized
                mR1_theory = has_m_Reporters1 ? df.m_Reporters1[1] : nothing
                mR2_theory = has_mReporters2 ? df.mReporters2[1] : nothing
            elseif has_cc_Reporters
                # cc_Reporters from write_correlation_functions is raw uncentered (E[xy]), use it directly
                cc_Reporters_theory_unnorm_all = df.cc_Reporters
                mR1_theory = has_m_Reporters1 ? df.m_Reporters1[1] : nothing
                mR2_theory = has_mReporters2 ? df.mReporters2[1] : nothing
            end

            # Read Reporter autocovariances if available
            if has_ac1_Reporters
                ac1_Reporters_theory_all = df.ac1_Reporters
            end
            if has_ac2_Reporters
                ac2_Reporters_theory_all = df.ac2_Reporters
            end

            # Interpolate theoretical values to match empirical lags
            function interpolate_theory(theory_all, tau_data, lags)
                theory_interp = Float64[]
                rtol = 1e-10  # Relative tolerance for float comparison
                for lag in lags
                    # Try to find exact match using approximate equality
                    idx = findfirst(x -> isapprox(x, lag, rtol=rtol), tau_data)
                    if !isnothing(idx)
                        push!(theory_interp, theory_all[idx])
                    else
                        # Linear interpolation
                        idx1 = findlast(x -> x < lag, tau_data)
                        idx2 = findfirst(x -> x > lag, tau_data)
                        if !isnothing(idx1) && !isnothing(idx2)
                            w1 = (tau_data[idx2] - lag) / (tau_data[idx2] - tau_data[idx1])
                            w2 = (lag - tau_data[idx1]) / (tau_data[idx2] - tau_data[idx1])
                            val = w1 * theory_all[idx1] + w2 * theory_all[idx2]
                            push!(theory_interp, val)
                        else
                            @warn "Could not find or interpolate lag $lag, skipping"
                            return nothing
                        end
                    end
                end
                return length(theory_interp) == length(lags) ? Vector{Float64}(theory_interp) : nothing
            end

            # Interpolate uncentered theory to match the filtered lags (same lags used for empirical)
            cc_ON_theory_unnorm = interpolate_theory(cc_ON_theory_unnorm_all, tau_data, lags)
            if isnothing(cc_ON_theory_unnorm)
                @warn "Could not interpolate ON state theory to match lags, skipping"
                continue
            end

            # Transform theory to match empirical CorrelationTrait
            # Theory is raw uncentered unnormalized: E[XY]
            # Apply the same transformations that were applied to empirical data
            if isnothing(mON1_theory) || isnothing(mON2_theory)
                @warn "Missing theoretical means (mON1_theory or mON2_theory), cannot transform theory. Skipping."
                continue
            end
            
            # Apply centering and normalization based on empirical CorrelationTrait
            needs_centering = empirical_correlation_algorithm.centering != :none
            needs_normalization = empirical_correlation_algorithm.normalization != :none
            
            # Compute normalization factor using global means (for global_mean normalization)
            # For windowed_mean or per_trace_mean normalization, this would need windowed/per-trace means
            # but for scoring we use global means as approximation
            norm_factor_ON = (mON1_theory * mON2_theory) > 0 ? (mON1_theory * mON2_theory) : 1.0
            
            if needs_centering
                cc_ON_theory_centered = cc_ON_theory_unnorm .- norm_factor_ON
            else
                cc_ON_theory_centered = cc_ON_theory_unnorm
            end
            
            if needs_normalization
                cc_ON_theory = cc_ON_theory_centered ./ norm_factor_ON
            else
                cc_ON_theory = cc_ON_theory_centered
            end

            # Empirical is already transformed according to CorrelationTrait
            cc_ON_empirical = Vector{Float64}(data_result.cc_unnormalized)

            # Verify all lengths match (they should by design)
            if length(cc_ON_theory) != length(lags)
                error("ON state theory length mismatch: theory=$(length(cc_ON_theory)), lags=$(length(lags))")
            end
            if length(cc_ON_empirical) != length(lags)
                error("ON state empirical length mismatch: empirical=$(length(cc_ON_empirical)), lags=$(length(lags))")
            end

            # Compute ON state scores (comparing normalized values)
            # Get standard error if available (convert to Vector{Float64} if not nothing)
            cc_ON_se = hasproperty(data_result, :cc_unnormalized_se) && !isnothing(data_result.cc_unnormalized_se) ? Vector{Float64}(data_result.cc_unnormalized_se) : nothing
            scores_ON = compute_score_metrics(Vector{Float64}(cc_ON_theory), Vector{Float64}(cc_ON_empirical), cc_ON_se)

            # Compute ON state autocovariance scores if available
            scores_ac1_ON = nothing
            scores_ac2_ON = nothing
            ac1_ON_theory = nothing
            ac1_ON_empirical = nothing
            ac2_ON_theory = nothing
            ac2_ON_empirical = nothing
            if !isnothing(ac1_ON_theory_all) && hasproperty(data_result, :ac1_unnormalized)
                ac1_ON_theory_unnorm = interpolate_theory(ac1_ON_theory_all, tau_data, lags)
                if !isnothing(ac1_ON_theory_unnorm) && length(ac1_ON_theory_unnorm) == length(lags)
                    # Transform theory using empirical CorrelationTrait
                    if isnothing(mON1_theory)
                        @warn "Missing mON1_theory, cannot transform ac1_ON theory. Skipping."
                    else
                        norm_factor_ac1 = mON1_theory^2 > 0 ? mON1_theory^2 : 1.0
                        if needs_centering
                            ac1_ON_theory_centered = ac1_ON_theory_unnorm .- norm_factor_ac1
                        else
                            ac1_ON_theory_centered = ac1_ON_theory_unnorm
                        end
                        if needs_normalization
                            ac1_ON_theory = ac1_ON_theory_centered ./ norm_factor_ac1
                        else
                            ac1_ON_theory = ac1_ON_theory_centered
                        end
                        ac1_ON_empirical = Vector{Float64}(data_result.ac1_unnormalized)
                        ac1_ON_se = hasproperty(data_result, :ac1_unnormalized_se) && !isnothing(data_result.ac1_unnormalized_se) ? Vector{Float64}(data_result.ac1_unnormalized_se) : nothing
                        scores_ac1_ON = compute_score_metrics(Vector{Float64}(ac1_ON_theory), ac1_ON_empirical, ac1_ON_se)
                    end
                end
            end
            if !isnothing(ac2_ON_theory_all) && hasproperty(data_result, :ac2_unnormalized)
                ac2_ON_theory_unnorm = interpolate_theory(ac2_ON_theory_all, tau_data, lags)
                if !isnothing(ac2_ON_theory_unnorm) && length(ac2_ON_theory_unnorm) == length(lags)
                    # Transform theory using empirical CorrelationTrait
                    if isnothing(mON2_theory)
                        @warn "Missing mON2_theory, cannot transform ac2_ON theory. Skipping."
                    else
                        norm_factor_ac2 = mON2_theory^2 > 0 ? mON2_theory^2 : 1.0
                        if needs_centering
                            ac2_ON_theory_centered = ac2_ON_theory_unnorm .- norm_factor_ac2
                        else
                            ac2_ON_theory_centered = ac2_ON_theory_unnorm
                        end
                        if needs_normalization
                            ac2_ON_theory = ac2_ON_theory_centered ./ norm_factor_ac2
                        else
                            ac2_ON_theory = ac2_ON_theory_centered
                        end
                        ac2_ON_empirical = Vector{Float64}(data_result.ac2_unnormalized)
                        ac2_ON_se = hasproperty(data_result, :ac2_unnormalized_se) && !isnothing(data_result.ac2_unnormalized_se) ? Vector{Float64}(data_result.ac2_unnormalized_se) : nothing
                        scores_ac2_ON = compute_score_metrics(Vector{Float64}(ac2_ON_theory), ac2_ON_empirical, ac2_ON_se)
                    end
                end
            end

            # Reporter scores (if available) - use same lags
            scores_Reporters = nothing
            cc_Reporters_theory = nothing
            cc_Reporters_empirical = nothing
            if !isnothing(cc_Reporters_theory_unnorm_all)
                cc_Reporters_theory_unnorm = interpolate_theory(cc_Reporters_theory_unnorm_all, tau_data, lags)
                if !isnothing(cc_Reporters_theory_unnorm)
                    # Transform theory using empirical CorrelationTrait
                    if isnothing(mR1_theory) || isnothing(mR2_theory)
                        @warn "Missing theoretical reporter means (mR1_theory or mR2_theory), cannot transform theory. Skipping Reporter scores."
                    else
                        norm_factor_Reporters = (mR1_theory * mR2_theory) > 0 ? (mR1_theory * mR2_theory) : 1.0
                        if needs_centering
                            cc_Reporters_theory_centered = cc_Reporters_theory_unnorm .- norm_factor_Reporters
                        else
                            cc_Reporters_theory_centered = cc_Reporters_theory_unnorm
                        end
                        if needs_normalization
                        cc_Reporters_theory = cc_Reporters_theory_centered ./ norm_factor_Reporters
                    else
                        cc_Reporters_theory = cc_Reporters_theory_centered
                    end

                        # Empirical is already centered and normalized
                        cc_Reporters_empirical = Vector{Float64}(reporter_result.cc)

                        # Verify all lengths match (they should by design)
                        if length(cc_Reporters_empirical) != length(lags)
                            @warn "Reporter empirical length mismatch: empirical=$(length(cc_Reporters_empirical)), lags=$(length(lags)). Skipping Reporter scores."
                        elseif length(cc_Reporters_theory) != length(lags)
                            @warn "Reporter theory length mismatch: theory=$(length(cc_Reporters_theory)), lags=$(length(lags)). Skipping Reporter scores."
                        else
                            # Compute Reporter scores (comparing normalized values)
                            # Get standard error if available (convert to Vector{Float64} if not nothing)
                            cc_Reporters_se = hasproperty(reporter_result, :cc_se) && !isnothing(reporter_result.cc_se) ? Vector{Float64}(reporter_result.cc_se) : nothing
                            scores_Reporters = compute_score_metrics(Vector{Float64}(cc_Reporters_theory), cc_Reporters_empirical, cc_Reporters_se)
                        end
                    end
                end
            end

            # Compute Reporter autocovariance scores if available
            scores_ac1_Reporters = nothing
            scores_ac2_Reporters = nothing
            ac1_Reporters_theory = nothing
            ac1_Reporters_empirical = nothing
            ac2_Reporters_theory = nothing
            ac2_Reporters_empirical = nothing
            if !isnothing(ac1_Reporters_theory_all) && hasproperty(reporter_result, :ac1)
                ac1_Reporters_theory_unnorm = interpolate_theory(ac1_Reporters_theory_all, tau_data, lags)
                if !isnothing(ac1_Reporters_theory_unnorm) && length(ac1_Reporters_theory_unnorm) == length(lags)
                    # Transform theory using empirical CorrelationTrait
                    if isnothing(mR1_theory)
                        @warn "Missing mR1_theory, cannot transform ac1_Reporters theory. Skipping."
                    else
                        norm_factor_ac1R = mR1_theory^2 > 0 ? mR1_theory^2 : 1.0
                        if needs_centering
                            ac1_Reporters_theory_centered = ac1_Reporters_theory_unnorm .- norm_factor_ac1R
                        else
                            ac1_Reporters_theory_centered = ac1_Reporters_theory_unnorm
                        end
                        if needs_normalization
                            ac1_Reporters_theory = ac1_Reporters_theory_centered ./ norm_factor_ac1R
                        else
                            ac1_Reporters_theory = ac1_Reporters_theory_centered
                        end
                        ac1_Reporters_empirical = Vector{Float64}(reporter_result.ac1)
                        ac1_Reporters_se = hasproperty(reporter_result, :ac1_se) && !isnothing(reporter_result.ac1_se) ? Vector{Float64}(reporter_result.ac1_se) : nothing
                        scores_ac1_Reporters = compute_score_metrics(Vector{Float64}(ac1_Reporters_theory), ac1_Reporters_empirical, ac1_Reporters_se)
                    end
                end
            end
            if !isnothing(ac2_Reporters_theory_all) && hasproperty(reporter_result, :ac2)
                ac2_Reporters_theory_unnorm = interpolate_theory(ac2_Reporters_theory_all, tau_data, lags)
                if !isnothing(ac2_Reporters_theory_unnorm) && length(ac2_Reporters_theory_unnorm) == length(lags)
                    # Transform theory using empirical CorrelationTrait
                    if isnothing(mR2_theory)
                        @warn "Missing mR2_theory, cannot transform ac2_Reporters theory. Skipping."
                    else
                        norm_factor_ac2R = mR2_theory^2 > 0 ? mR2_theory^2 : 1.0
                        if needs_centering
                            ac2_Reporters_theory_centered = ac2_Reporters_theory_unnorm .- norm_factor_ac2R
                        else
                            ac2_Reporters_theory_centered = ac2_Reporters_theory_unnorm
                        end
                        if needs_normalization
                            ac2_Reporters_theory = ac2_Reporters_theory_centered ./ norm_factor_ac2R
                        else
                            ac2_Reporters_theory = ac2_Reporters_theory_centered
                        end
                        ac2_Reporters_empirical = Vector{Float64}(reporter_result.ac2)
                        ac2_Reporters_se = hasproperty(reporter_result, :ac2_se) && !isnothing(reporter_result.ac2_se) ? Vector{Float64}(reporter_result.ac2_se) : nothing
                        scores_ac2_Reporters = compute_score_metrics(Vector{Float64}(ac2_Reporters_theory), ac2_Reporters_empirical, ac2_Reporters_se)
                    end
                end
            end

            # Store results (all unnormalized)
            result_dict = Dict(
                :coupling_model => coupling_model_str,
                :model => model_str,
                :cc_ON_empirical => cc_ON_empirical,
                :cc_ON_theory => cc_ON_theory,
                :cc_ON_l2_norm => scores_ON.l2_norm,
                :cc_ON_linf_norm => scores_ON.linf_norm,
                :mON1_empirical => data_result.mean1,
                :mON2_empirical => data_result.mean2,
                :mON1_theory => mON1_theory,
                :mON2_theory => mON2_theory,
                :lags => lags
            )

            # Add weighted scores and chi-squared if available
            if hasproperty(scores_ON, :weighted_l2_norm)
                result_dict[:cc_ON_weighted_l2_norm] = scores_ON.weighted_l2_norm
                result_dict[:cc_ON_weighted_linf_norm] = scores_ON.weighted_linf_norm
                result_dict[:cc_ON_chi_squared] = scores_ON.chi_squared
                result_dict[:cc_ON_reduced_chi_squared] = scores_ON.reduced_chi_squared
            end

            # Add ON state autocovariance results if available
            if !isnothing(scores_ac1_ON)
                result_dict[:ac1_ON_empirical] = ac1_ON_empirical
                result_dict[:ac1_ON_theory] = ac1_ON_theory
                result_dict[:ac1_ON_l2_norm] = scores_ac1_ON.l2_norm
                result_dict[:ac1_ON_linf_norm] = scores_ac1_ON.linf_norm
                if hasproperty(scores_ac1_ON, :weighted_l2_norm)
                    result_dict[:ac1_ON_weighted_l2_norm] = scores_ac1_ON.weighted_l2_norm
                    result_dict[:ac1_ON_weighted_linf_norm] = scores_ac1_ON.weighted_linf_norm
                    result_dict[:ac1_ON_chi_squared] = scores_ac1_ON.chi_squared
                    result_dict[:ac1_ON_reduced_chi_squared] = scores_ac1_ON.reduced_chi_squared
                end
            end
            if !isnothing(scores_ac2_ON)
                result_dict[:ac2_ON_empirical] = ac2_ON_empirical
                result_dict[:ac2_ON_theory] = ac2_ON_theory
                result_dict[:ac2_ON_l2_norm] = scores_ac2_ON.l2_norm
                result_dict[:ac2_ON_linf_norm] = scores_ac2_ON.linf_norm
                if hasproperty(scores_ac2_ON, :weighted_l2_norm)
                    result_dict[:ac2_ON_weighted_l2_norm] = scores_ac2_ON.weighted_l2_norm
                    result_dict[:ac2_ON_weighted_linf_norm] = scores_ac2_ON.weighted_linf_norm
                    result_dict[:ac2_ON_chi_squared] = scores_ac2_ON.chi_squared
                    result_dict[:ac2_ON_reduced_chi_squared] = scores_ac2_ON.reduced_chi_squared
                end
            end

            # Add Reporter results if available (all unnormalized)
            if !isnothing(scores_Reporters) && !isnothing(cc_Reporters_theory) && !isnothing(cc_Reporters_empirical)
                result_dict[:cc_Reporters_empirical] = cc_Reporters_empirical
                result_dict[:cc_Reporters_theory] = cc_Reporters_theory
                result_dict[:cc_Reporters_l2_norm] = scores_Reporters.l2_norm
                result_dict[:cc_Reporters_linf_norm] = scores_Reporters.linf_norm
                if hasproperty(scores_Reporters, :weighted_l2_norm)
                    result_dict[:cc_Reporters_weighted_l2_norm] = scores_Reporters.weighted_l2_norm
                    result_dict[:cc_Reporters_weighted_linf_norm] = scores_Reporters.weighted_linf_norm
                    result_dict[:cc_Reporters_chi_squared] = scores_Reporters.chi_squared
                    result_dict[:cc_Reporters_reduced_chi_squared] = scores_Reporters.reduced_chi_squared
                end
                result_dict[:mR1_empirical] = reporter_result.mean1
                result_dict[:mR2_empirical] = reporter_result.mean2
                result_dict[:mR1_theory] = mR1_theory
                result_dict[:mR2_theory] = mR2_theory
            end

            # Add Reporter autocovariance results if available
            if !isnothing(scores_ac1_Reporters)
                result_dict[:ac1_Reporters_empirical] = ac1_Reporters_empirical
                result_dict[:ac1_Reporters_theory] = ac1_Reporters_theory
                result_dict[:ac1_Reporters_l2_norm] = scores_ac1_Reporters.l2_norm
                result_dict[:ac1_Reporters_linf_norm] = scores_ac1_Reporters.linf_norm
                if hasproperty(scores_ac1_Reporters, :weighted_l2_norm)
                    result_dict[:ac1_Reporters_weighted_l2_norm] = scores_ac1_Reporters.weighted_l2_norm
                    result_dict[:ac1_Reporters_weighted_linf_norm] = scores_ac1_Reporters.weighted_linf_norm
                    result_dict[:ac1_Reporters_chi_squared] = scores_ac1_Reporters.chi_squared
                    result_dict[:ac1_Reporters_reduced_chi_squared] = scores_ac1_Reporters.reduced_chi_squared
                end
            end
            if !isnothing(scores_ac2_Reporters)
                result_dict[:ac2_Reporters_empirical] = ac2_Reporters_empirical
                result_dict[:ac2_Reporters_theory] = ac2_Reporters_theory
                result_dict[:ac2_Reporters_l2_norm] = scores_ac2_Reporters.l2_norm
                result_dict[:ac2_Reporters_linf_norm] = scores_ac2_Reporters.linf_norm
                if hasproperty(scores_ac2_Reporters, :weighted_l2_norm)
                    result_dict[:ac2_Reporters_weighted_l2_norm] = scores_ac2_Reporters.weighted_l2_norm
                    result_dict[:ac2_Reporters_weighted_linf_norm] = scores_ac2_Reporters.weighted_linf_norm
                    result_dict[:ac2_Reporters_chi_squared] = scores_ac2_Reporters.chi_squared
                    result_dict[:ac2_Reporters_reduced_chi_squared] = scores_ac2_Reporters.reduced_chi_squared
                end
            end

            results[results_key] = NamedTuple(result_dict)

        catch e
            @warn "Error reading crosscovariance file '$crosscov_file': $e, skipping"
            continue
        end
    end

    return results
end

"""
    compute_score_metrics(theory::Vector{Float64}, empirical::Vector{Float64}, se::Union{Vector{Float64}, Nothing}=nothing)

Compute scoring metrics comparing theoretical and empirical cross-covariances.

# Arguments
- `theory::Vector{Float64}`: Theoretical cross-covariance values
- `empirical::Vector{Float64}`: Empirical cross-covariance values  
- `se::Union{Vector{Float64}, Nothing}`: Standard error of empirical values (optional)

# Returns
- NamedTuple with:
  - `l2_norm`: Standard L² norm (unweighted)
  - `linf_norm`: Standard L∞ norm (unweighted)
  - `weighted_l2_norm`: Weighted L² norm (if SE provided)
  - `weighted_linf_norm`: Weighted L∞ norm (max z-score, if SE provided)
  - `chi_squared`: Chi-squared statistic (if SE provided)
  - `reduced_chi_squared`: Reduced chi-squared (chi²/dof, if SE provided)
"""
function compute_score_metrics(theory::Vector{Float64}, empirical::Vector{Float64}, se::Union{Vector{Float64},Nothing}=nothing)
    diff = theory .- empirical

    # Standard L² and L∞ norms (unweighted)
    l2_norm = sqrt(sum(diff .^ 2) / length(diff))
    linf_norm = maximum(abs.(diff))

    result = (l2_norm=l2_norm, linf_norm=linf_norm)

    # If uncertainty (standard error) is provided, compute weighted scores
    if !isnothing(se) && length(se) == length(diff)
        # Avoid division by zero: use a minimum SE threshold
        se_safe = max.(se, 1e-10)

        # Weighted L² norm: sqrt(sum((diff ./ se_safe) .^ 2) / length(diff))
        # This is like the square root of the mean chi-squared per data point
        weighted_l2_norm = sqrt(sum((diff ./ se_safe) .^ 2) / length(diff))

        # Weighted L∞ norm: max(|theory - empirical| / se)
        # This is the maximum z-score
        weighted_linf_norm = maximum(abs.(diff ./ se_safe))

        # Chi-squared statistic: sum((diff ./ se_safe) .^ 2)
        chi_squared = sum((diff ./ se_safe) .^ 2)

        # Reduced chi-squared (chi-squared per degree of freedom)
        reduced_chi_squared = chi_squared / length(diff)

        return merge(result, (
            weighted_l2_norm=weighted_l2_norm,
            weighted_linf_norm=weighted_linf_norm,
            chi_squared=chi_squared,
            reduced_chi_squared=reduced_chi_squared
        ))
    end

    return result
end

"""
    summarize_model_scores(results::Dict{String, NamedTuple})

Print a summary table of model scores sorted by L² norm (best fit first).
"""
function summarize_model_scores(results::Dict{String,NamedTuple}; sort_by::Symbol=:cc_ON_l2_norm, use_reporters::Bool=false, csv_file::Union{String,Nothing}=nothing)
    # Check what types of scores are available
    first_result = first(values(results))
    has_ON = hasproperty(first_result, :cc_ON_l2_norm)
    has_Reporters = hasproperty(first_result, :cc_Reporters_l2_norm)
    has_ON_weighted = hasproperty(first_result, :cc_ON_weighted_l2_norm)
    has_Reporters_weighted = hasproperty(first_result, :cc_Reporters_weighted_l2_norm)

    # Determine sort key
    sort_key = if sort_by == :cc_ON_l2_norm && has_ON
        k -> get(results[k], :cc_ON_l2_norm, Inf)
    elseif sort_by == :cc_Reporters_l2_norm && has_Reporters
        k -> get(results[k], :cc_Reporters_l2_norm, Inf)
    elseif sort_by == :cc_ON_weighted_l2_norm && has_ON_weighted
        k -> get(results[k], :cc_ON_weighted_l2_norm, Inf)
    elseif sort_by == :cc_Reporters_weighted_l2_norm && has_Reporters_weighted
        k -> get(results[k], :cc_Reporters_weighted_l2_norm, Inf)
    elseif use_reporters && has_Reporters
        k -> get(results[k], :cc_Reporters_l2_norm, Inf)
    elseif has_ON
        k -> get(results[k], :cc_ON_l2_norm, Inf)
    else
        k -> Inf
    end

    println("\nModel Scoring Summary (sorted by $sort_by, best fit first):")
    println("="^100)

    # Build header
    header_cols = ["Model"]
    if has_ON
        append!(header_cols, ["ON_L2", "ON_Linf"])
        if has_ON_weighted
            append!(header_cols, ["ON_wL2", "ON_wLinf", "ON_chi2", "ON_chi2_dof"])
        end
        append!(header_cols, ["ON_m1", "ON_m2"])
    end
    if has_Reporters
        append!(header_cols, ["Rep_L2", "Rep_Linf"])
        if has_Reporters_weighted
            append!(header_cols, ["Rep_wL2", "Rep_wLinf", "Rep_chi2", "Rep_chi2_dof"])
        end
        append!(header_cols, ["Rep_m1", "Rep_m2"])
    end

    header_str = @sprintf("%-10s", header_cols[1])
    for i in 2:length(header_cols)
        header_str *= @sprintf(" %12s", header_cols[i])
    end
    println(header_str)
    println("-"^100)

    sorted_models = sort(collect(keys(results)), by=sort_key)

    # Prepare CSV data if requested
    csv_data = Dict{String,Vector{Any}}()
    if !isnothing(csv_file)
        csv_data["Model"] = String[]
        for col in header_cols[2:end]
            csv_data[col] = Float64[]
        end
    end

    for model in sorted_models
        r = results[model]
        row_values = Any[model]

        if has_ON
            push!(row_values, get(r, :cc_ON_l2_norm, NaN))
            push!(row_values, get(r, :cc_ON_linf_norm, NaN))
            if has_ON_weighted
                push!(row_values, get(r, :cc_ON_weighted_l2_norm, NaN))
                push!(row_values, get(r, :cc_ON_weighted_linf_norm, NaN))
                push!(row_values, get(r, :cc_ON_chi_squared, NaN))
                push!(row_values, get(r, :cc_ON_reduced_chi_squared, NaN))
            end
            push!(row_values, get(r, :mON1_empirical, NaN))
            push!(row_values, get(r, :mON2_empirical, NaN))
        end

        if has_Reporters
            push!(row_values, get(r, :cc_Reporters_l2_norm, NaN))
            push!(row_values, get(r, :cc_Reporters_linf_norm, NaN))
            if has_Reporters_weighted
                push!(row_values, get(r, :cc_Reporters_weighted_l2_norm, NaN))
                push!(row_values, get(r, :cc_Reporters_weighted_linf_norm, NaN))
                push!(row_values, get(r, :cc_Reporters_chi_squared, NaN))
                push!(row_values, get(r, :cc_Reporters_reduced_chi_squared, NaN))
            end
            push!(row_values, get(r, :mR1_empirical, NaN))
            push!(row_values, get(r, :mR2_empirical, NaN))
        end

        # Print formatted row
        row_str = @sprintf("%-10s", row_values[1])
        for i in 2:length(row_values)
            row_str *= @sprintf(" %12.6f", row_values[i])
        end
        println(row_str)

        # Populate CSV data
        if !isnothing(csv_file)
            push!(csv_data["Model"], string(row_values[1]))
            for (i, col) in enumerate(header_cols[2:end])
                push!(csv_data[col], Float64(row_values[i+1]))
            end
        end
    end

    println("="^100)
    println("\nNote: Empirical values are the same for all models (computed from data).")
    println("      Differences are in theoretical predictions, which determine the norms.")
    println()

    # Write to CSV if requested
    if !isnothing(csv_file)
        df = DataFrame(csv_data)
        CSV.write(csv_file, df)
        println("Scoring summary saved to '$csv_file'")
    end

    return nothing
end




# ============================================================================
# PLOTTING FUNCTIONS (OBSOLETE - Plots dependency removed)
# ============================================================================
# Functions: plot_traces, plot_histogram (multiple overloads)
# Creation date: OLD (pre-2025, user-created)
# Status: Kept for assessment, but Plots.jl dependency removed for Biowulf compatibility
# TODO: Remove or move to separate visualization module

# """
#     plot_traces(df::DataFrame; 
#                cols::Vector{Int}=collect(1:3),  # Convert UnitRange to Vector
#                with_errors::Bool=true,
#                normalize::Bool=true,
#                title::String="Traces",
#                alpha::Float64=0.3)

# Plot multiple traces in separate panels for data and model.

# # Arguments
# - `df::DataFrame`: DataFrame containing data/model columns
# - `cols::Vector{Int}`: Column numbers to plot (default: [1,2,3])
# - `with_errors::Bool`: Include error ribbons (default: true)
# - `normalize::Bool`: Normalize data (default: true)
# - `title::String`: Plot title (default: "Traces")
# - `alpha::Float64`: Transparency for ribbons (default: 0.3)

# # Returns
# - `Plots.Plot`: Combined plot object
# """
# function plot_traces(df::DataFrame;
#     cols::Vector{Int}=collect(1:2),  # Changed this line
#     with_errors::Bool=false,
#     normalize::Bool=false,
#     title::String="Traces",
#     alpha::Float64=0.3)

#     n = length(cols)
#     p = plot(layout=(n, 1), legend=:topleft, title=title)

#     for (j, i) in enumerate(cols)
#         data_col = Symbol("data$i")
#         mean_col = Symbol("model_mean$i")
#         std_col = Symbol("model_std$i")

#         # Get data, filtering out missing values
#         y_data = df[!, data_col]
#         y_model = df[!, mean_col]
#         y_std = hasproperty(df, std_col) ? df[!, std_col] : nothing

#         # Filter out missing values
#         valid_indices = .!ismissing.(y_data) .& .!ismissing.(y_model)
#         y_data = y_data[valid_indices]
#         y_model = y_model[valid_indices]
#         y_std = y_std !== nothing ? y_std[valid_indices] : nothing

#         # Normalize if needed
#         if normalize
#             y_data ./= maximum(y_data)
#             y_model ./= maximum(y_model)
#             y_std = y_std !== nothing ? y_std ./ maximum(y_model) : nothing
#         end

#         # Plot data points in the first row
#         scatter!(p[j, 1], y_data, label="Data $i", markersize=3)

#         # Plot model with optional error ribbon in the second row
#         if with_errors && y_std !== nothing
#             plot!(p[j, 1], y_model, err=y_std, label="Model $i", alpha=alpha)
#         else
#             plot!(p[j, 1], y_model, label="Model $i")
#         end
#     end

#     # xlabel!(p[1, :], "Time point")
#     # ylabel!(p[1, :], normalize ? "Normalized value" : "Value")
#     # xlabel!(p[2, :], "Time point")
#     # ylabel!(p[2, :], normalize ? "Normalized value" : "Value")

#     return p
# end


# """
#     plot_histogram(ratefile::String, datapath; root=".", row=2)

#     plot_histogram()
#     plot_histogram(ratefile::String,datapath;fish=false,root=".",row=2)

#     functions to plot data and model predicted histograms

# """
# function plot_histogram(ratefile::String, datapath; root=".", row=2)
#     fish = false
#     r = readrow(ratefile, row)
#     println("r: ", r)
#     parts = fields(ratefile)
#     label = parts.label
#     cond = parts.cond
#     G = parts.model
#     data = data_rna(parts.gene, parts.cond, datapath, fish, parts.label, root)
#     model = model_rna(r, r, parse(Int, parts.model), parse(Int, parts.nalleles), 0.01, [], (), 0)
#     # model_rna(r,[rateprior[i]],G,nalleles,cv,[fittedparam[i]],fixedeffects,0)
#     plot_histogram(data, model)
#     return data, model
# end

# function plot_histogram(gene::String, cell::String, G::Int, cond::String, ratefile::String, datapath::String, root::String=".")
#     fish = false
#     rates = readdlm(ratefile, ',', header=true)
#     r = rates[1][findfirst(rates[1][:, 1] .== gene)[1], 2:end]
#     data = data_rna(gene, cond, datapath, fish, "label", root)
#     nalleles = alleles(gene, cell, root)
#     model = model_rna(r, [], G, nalleles, 0.01, [], (), 0)
#     println("typeof(model): ", typeof(model))
#     println("typeof(data): ", typeof(data))
#     m = plot_histogram(data, model)
#     return m, data, model
# end

# function plot_histogram(gene::String, cell::String, G::String, cond::String, label::String, ratefolder::String, datapath::String, nsets::Int, root::String, fittedparam=[1], verbose=false)
#     fish = false
#     data = data_rna(gene, cond, datapath, fish, label, root)
#     model = model_rna(gene, cell, G, fish, 0.01, fittedparam, (), label, ratefolder, nsets, root, data, verbose)
#     m = plot_histogram(data, model)
#     return m, data, model
# end

# function plot_histogram(data::AbstractRNAData{Array{Array,1}}, model)
#     h = predictedarray(model.rates, data, model)
#     figure(data.gene)
#     for i in eachindex(h)
#         plot(h[i])
#         plot(normalize_histogram(data.histRNA[i]))
#         savefig(string(i))
#     end
#     return h
# end
# function plot_histogram(data::AbstractRNAData{Array{Float64,1}}, model)
#     h = predictedfn(get_param(model), data, model)
#     plt = plot(h)
#     plot!(plt, normalize_histogram(data.histRNA))
#     display(plt)
#     return h
# end

# function plot_histogram(data::RNAOnOffData, model::AbstractGeneTransitionModel, filename="")
#     h = predictedarray(model.rates, data, model)
#     plt1 = plot(data.bins, h[1])
#     plot!(plt1, data.bins, normalize_histogram(data.OFF))
#     plt2 = plot(data.bins, h[2])
#     plot!(plt2, data.bins, normalize_histogram(data.ON))
#     plt3 = plot(h[3])
#     plot!(plt3, normalize_histogram(data.histRNA))
#     plt = plot(plt1, plt2, plt3, layout=(3, 1))
#     display(plt)
#     if ~isempty(filename)
#         savefig(filename)
#     end
#     return h
# end

# function plot_histogram(data::TraceRNAData, model::AbstractGeneTransitionModel, filename="")
#     M = make_mat_M(model.components.mcomponents, model.rates)
#     h = steady_state(M, model.components.mcomponents.nT, model.nalleles, data.nRNA)
#     plt = plot(h)
#     plot!(plt, normalize_histogram(data.histRNA))
#     display(plt)
#     if ~isempty(filename)
#         savefig(filename)
#     end
#     return h
# end

# function plot_histogram(data::RNAData{T1,T2}, model::AbstractGMmodel, save=false) where {T1<:Array,T2<:Array}
#     m = predictedarray(model.rates, data, model)
#     for i in eachindex(m)
#         plt = plot(m[i])
#         plot!(normalize_histogram(data.histRNA[i]), show=true)
#         if save
#             savefig()
#         else
#             display(plt)
#         end
#     end
#     println("deviance: ", deviance(data, model))
#     println("loglikelihood: ", loglikelihood(get_param(model), data, model)[1])
#     return m
# end

# function plot_histogram(data::RNAData, model::AbstractGMmodel)
#     h = predictedfn(get_param(model), data, model)
#     plt = plot(h)
#     plot!(normalize_histogram(data.histRNA))
#     display(plt)
#     return h
# end
