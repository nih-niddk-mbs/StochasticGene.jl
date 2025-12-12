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
    for (root, dirs, files) in walkdir(folder)
        for f in files
            if occursin("rates", f) 
                write_trace_dataframe(joinpath(root, f), datapath, interval, ratetype, start, stop, probfn, noiseparams, splicetype, hlabel=hlabel, state=state, grid=grid, zeromedian=zeromedian, datacol=datacol)
            end
        end
    end
end



"""
    write_trace_dataframe(file, datapath::String, datacond, interval, ratetype::String="median", start=1, stop=-1, probfn=prob_Gaussian, noiseparams=4, splicetype=""; hlabel="-h", state=true, coupling=tuple())

TBW
"""
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
                ac1, ac2, cc, ccON, tau, m1, m2, v1, v2, mON1, mON2, ac1ON, ac2ON = write_cov_file(file, transitions, G, R, S, insertstep, interval, pattern, lags, probfn, ratetype)
                out = replace(file, "rates" => "crosscovariance", ".txt" => ".csv")
                ac1 = ac1/v1
                ac2 = ac2/v2
                ac1ON = ac1ON/mON1^2
                ac2ON = ac2ON/mON2^2
                CSV.write(out, DataFrame(tau=tau, crosscorrelation=cc/sqrt(v1*v2), autocor1=[reverse(ac1); ac1[2:end]], autocor2=[reverse(ac2); ac2[2:end]], cc_ON=ccON/mON1/mON2, ac1_ON=[reverse(ac1ON); ac1ON[2:end]], ac2_ON=[reverse(ac2ON); ac2ON[2:end]]))
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
- Named tuple with fields:
    - `linf_norm`: L-infinity norms between theory and simulation (one per trial)
    - `l2_norm`: L2 norms between theory and simulation (one per trial)
    - `cc_theory`: Theoretical cross-covariance function
    - `cc_mean`: Mean cross-covariance from simulations
    - `cc_se`: Standard error of cross-covariance mean
    - `lags`: Time lags
    - `ac1_mean`: Mean autocovariance for first trace from simulations
    - `ac1_theory`: Theoretical autocovariance for first trace
    - `ac2_mean`: Mean autocovariance for second trace from simulations
    - `ac2_theory`: Theoretical autocovariance for second trace

# Example
```julia
ntrials = 100
result = simulate_trials(r, trans, G, R, S, insertstep, coupling, ntrials)
plot(1:ntrials, result.l2_norm, ylabel="L2 Norm", xlabel="Number of Trials")
plot(result.lags, result.cc_mean, label="Simulated")
plot!(result.lags, result.cc_theory, label="Theoretical")
```
"""

function simulate_trials(r::Vector, transitions::Tuple, G, R, S, insertstep, coupling, ntrials, trial_time=720.0, lag=60, stride=1, probfn=prob_Gaussian)
    positive_lags = collect(0:stride:lag)
    # Get both intensity and ON/OFF cross-covariances and autocovariances
    ac1, ac2, cc, ccON, full_lags, m1, m2, v1, v2, mON1, mON2, ac1ON, ac2ON = StochasticGene.covariance_functions(r, transitions, G, R, S, insertstep, 1.0, probfn, coupling, positive_lags)
    # Use ccON (ON/OFF cross-covariance) and ac1ON/ac2ON (ON/OFF autocovariances) for comparison
    simulate_trials(ac1, ac2, cc, ac1ON, ac2ON, ccON, mON1, mON2, v1, v2, r, transitions, G, R, S, insertstep, coupling, positive_lags, full_lags, ntrials, trial_time)
end

function simulate_trials(ac1, ac2, cc, ac1ON, ac2ON, ccON, mON1, mON2, v1, v2, r::Vector, transitions::Tuple, G, R, S, insertstep, coupling, lags_ac, lags, ntrials=1, trial_time::Float64=720.0)
    # Store per-trial results for parallel processing and bootstrap
    # Need to compute 6 covariance sets: intensity (AC1, AC2, CC) and ON state (AC1, AC2, CC)
    mean1_trials = Vector{Float64}(undef, ntrials)
    mean2_trials = Vector{Float64}(undef, ntrials)
    mean1_intensity_trials = Vector{Float64}(undef, ntrials)  # Intensity means
    mean2_intensity_trials = Vector{Float64}(undef, ntrials)  # Intensity means
    cc_raw_trials = Vector{Vector{Float64}}(undef, ntrials)  # ON state CC
    cc_intensity_raw_trials = Vector{Vector{Float64}}(undef, ntrials)  # Intensity CC
    ac1_raw_trials = Vector{Vector{Float64}}(undef, ntrials)  # Intensity AC
    ac2_raw_trials = Vector{Vector{Float64}}(undef, ntrials)  # Intensity AC
    ac1ON_raw_trials = Vector{Vector{Float64}}(undef, ntrials)  # ON state AC
    ac2ON_raw_trials = Vector{Vector{Float64}}(undef, ntrials)  # ON state AC
    trace1_trials = Vector{Vector{Float64}}(undef, ntrials)  # Store per-trial traces for variance
    trace2_trials = Vector{Vector{Float64}}(undef, ntrials)
    trace1_intensity_trials = Vector{Vector{Float64}}(undef, ntrials)  # Store intensity traces for AC
    trace2_intensity_trials = Vector{Vector{Float64}}(undef, ntrials)
    
    # Parallel processing: run trials concurrently
    # With 1 thread, this runs sequentially with minimal overhead
    Threads.@threads for n in 1:ntrials
        # Get traces with col=[2,3] returns matrix with columns: [intensity1, reporters1, intensity2, reporters2]
        t = simulate_trace_vector(r, transitions, G, R, S, insertstep, coupling, 1.0, trial_time+100., 1, col=[2,3])
        
        # Extract from the combined trace matrix (skip warmup period)
        # t[1] is a matrix with columns: [intensity1, reporters1, intensity2, reporters2]
        trace_data = t[1][100:end, :]
        
        # Extract intensity traces (columns 1 and 3)
        trace1_intensity = trace_data[:, 1]  # Intensity for component 1 (enhancer)
        trace2_intensity = trace_data[:, 3]  # Intensity for component 2 (gene)
        
        # Extract reporter traces (columns 2 and 4) and convert to binary ON/OFF
        reporters1 = trace_data[:, 2]  # Reporters for component 1 (enhancer)
        reporters2 = trace_data[:, 4]  # Reporters for component 2 (gene)
        trace1_ON = Float64.(reporters1 .> 0.0)  # Binary ON/OFF for enhancer (reporter > 0)
        trace2_ON = Float64.(reporters2 .> 0.0)  # Binary ON/OFF for gene (reporter > 0)
        
        # Store traces for variance calculation (no locking needed - each thread writes to different index)
        trace1_trials[n] = trace1_ON
        trace2_trials[n] = trace2_ON
        trace1_intensity_trials[n] = trace1_intensity
        trace2_intensity_trials[n] = trace2_intensity
        
        # Compute per-trial means for bootstrap
        mean1_trials[n] = mean(trace1_ON)  # ON state mean (fraction ON)
        mean2_trials[n] = mean(trace2_ON)  # ON state mean (fraction ON)
        mean1_intensity_trials[n] = mean(trace1_intensity)  # Intensity mean
        mean2_intensity_trials[n] = mean(trace2_intensity)  # Intensity mean
        
        # Compute raw cross-correlations and autocorrelations WITHOUT demeaning
        # ON state correlations
        cc_raw_trials[n] = StatsBase.crosscov(trace1_ON, trace2_ON, lags, demean=false)  # ON state CC
        ac1ON_raw_trials[n] = StatsBase.autocov(trace1_ON, lags_ac, demean=false)  # ON state AC1
        ac2ON_raw_trials[n] = StatsBase.autocov(trace2_ON, lags_ac, demean=false)  # ON state AC2
        # Intensity correlations
        cc_intensity_raw_trials[n] = StatsBase.crosscov(trace1_intensity, trace2_intensity, lags, demean=false)  # Intensity CC
        ac1_raw_trials[n] = StatsBase.autocov(trace1_intensity, lags_ac, demean=false)  # Intensity AC1
        ac2_raw_trials[n] = StatsBase.autocov(trace2_intensity, lags_ac, demean=false)  # Intensity AC2
    end
    
    # Compute all statistics from stored results (no online updates needed)
    # Combine all traces for overall mean and variance
    # ON state traces
    trace1_all = vcat(trace1_trials...)
    trace2_all = vcat(trace2_trials...)
    mean1 = mean(trace1_all)  # ON state mean (fraction ON)
    mean2 = mean(trace2_all)  # ON state mean (fraction ON)
    var1 = var(trace1_all)
    var2 = var(trace2_all)
    
    # Intensity traces
    trace1_intensity_all = vcat(trace1_intensity_trials...)
    trace2_intensity_all = vcat(trace2_intensity_trials...)
    mean1_intensity = mean(trace1_intensity_all)  # Intensity mean
    mean2_intensity = mean(trace2_intensity_all)  # Intensity mean
    var1_intensity = var(trace1_intensity_all)
    var2_intensity = var(trace2_intensity_all)
    
    # Compute mean of raw moments across trials (element-wise)
    # ON state correlations
    cc_raw_mean = [mean([cc_raw_trials[i][j] for i in 1:ntrials]) for j in 1:length(lags)]  # ON state CC
    ac1ON_raw_mean = [mean([ac1ON_raw_trials[i][j] for i in 1:ntrials]) for j in 1:length(lags_ac)]  # ON state AC1
    ac2ON_raw_mean = [mean([ac2ON_raw_trials[i][j] for i in 1:ntrials]) for j in 1:length(lags_ac)]  # ON state AC2
    # Intensity correlations
    cc_intensity_raw_mean = [mean([cc_intensity_raw_trials[i][j] for i in 1:ntrials]) for j in 1:length(lags)]  # Intensity CC
    ac1_raw_mean = [mean([ac1_raw_trials[i][j] for i in 1:ntrials]) for j in 1:length(lags_ac)]  # Intensity AC1
    ac2_raw_mean = [mean([ac2_raw_trials[i][j] for i in 1:ntrials]) for j in 1:length(lags_ac)]  # Intensity AC2
    
    # Compute variance (standard error) of raw moments across trials
    cc_var = [var([cc_raw_trials[i][j] for i in 1:ntrials]) for j in 1:length(lags)]  # ON state CC variance
    cc_intensity_var = [var([cc_intensity_raw_trials[i][j] for i in 1:ntrials]) for j in 1:length(lags)]  # Intensity CC variance
    ac1_var = [var([ac1_raw_trials[i][j] for i in 1:ntrials]) for j in 1:length(lags_ac)]  # Intensity AC1 variance
    ac2_var = [var([ac2_raw_trials[i][j] for i in 1:ntrials]) for j in 1:length(lags_ac)]  # Intensity AC2 variance
    ac1ON_var = [var([ac1ON_raw_trials[i][j] for i in 1:ntrials]) for j in 1:length(lags_ac)]  # ON state AC1 variance
    ac2ON_var = [var([ac2ON_raw_trials[i][j] for i in 1:ntrials]) for j in 1:length(lags_ac)]  # ON state AC2 variance
    
    # Standard error of the mean
    mean1_se = sqrt(var1 / length(trace1_all))
    mean2_se = sqrt(var2 / length(trace2_all))
    
    # Bootstrap confidence intervals for means (handles non-normal distributions)
    n_bootstrap = 10000
    # Use thread-local storage to avoid contention
    nthreads_mean = Threads.nthreads()
    mean1_bootstrap_threads = [Float64[] for _ in 1:nthreads_mean]
    mean2_bootstrap_threads = [Float64[] for _ in 1:nthreads_mean]
    Threads.@threads for _ in 1:n_bootstrap
        thread_id = Threads.threadid()
        # Resample with replacement from per-trial means
        bootstrap_sample = StatsBase.sample(mean1_trials, length(mean1_trials), replace=true)
        push!(mean1_bootstrap_threads[thread_id], mean(bootstrap_sample))
        bootstrap_sample = StatsBase.sample(mean2_trials, length(mean2_trials), replace=true)
        push!(mean2_bootstrap_threads[thread_id], mean(bootstrap_sample))
    end
    # Combine thread-local results
    mean1_bootstrap = vcat(mean1_bootstrap_threads...)
    mean2_bootstrap = vcat(mean2_bootstrap_threads...)
    mean1_lower = quantile(mean1_bootstrap, 0.025)
    mean1_median = median(mean1_bootstrap)
    mean1_upper = quantile(mean1_bootstrap, 0.975)
    mean2_lower = quantile(mean2_bootstrap, 0.025)
    mean2_median = median(mean2_bootstrap)
    mean2_upper = quantile(mean2_bootstrap, 0.975)
    
    # Keep SE for backward compatibility, but bootstrap intervals are more reliable for non-normal data
    mean1_se = length(trace1_all) > 1 ? sqrt(var1 / length(trace1_all)) : 0.0
    mean2_se = length(trace2_all) > 1 ? sqrt(var2 / length(trace2_all)) : 0.0
    
    # Subtract overall empirical means at the end to get cross-covariance and autocovariance
    # This ensures consistency: E[xy] - E[x]E[y] where all expectations are empirical
    # ON state correlations
    cc_mean = cc_raw_mean .- (mean1 * mean2)  # ON state CC
    ac1ON_mean = ac1ON_raw_mean .- (mean1^2)  # ON state AC1
    ac2ON_mean = ac2ON_raw_mean .- (mean2^2)  # ON state AC2
    # Intensity correlations
    cc_intensity_mean = cc_intensity_raw_mean .- (mean1_intensity * mean2_intensity)  # Intensity CC
    ac1_mean = ac1_raw_mean .- (mean1_intensity^2)  # Intensity AC1
    ac2_mean = ac2_raw_mean .- (mean2_intensity^2)  # Intensity AC2
    
    # Bootstrap confidence intervals for covariance functions (handles non-normal distributions)
    # For each lag, bootstrap resample trials and compute percentiles
    n_bootstrap_cov = 10000
    # Use thread-local storage to avoid contention
    nthreads = Threads.nthreads()
    cc_bootstrap_threads = [[Float64[] for _ in 1:length(lags)] for _ in 1:nthreads]  # ON state CC
    cc_intensity_bootstrap_threads = [[Float64[] for _ in 1:length(lags)] for _ in 1:nthreads]  # Intensity CC
    ac1_bootstrap_threads = [[Float64[] for _ in 1:length(lags_ac)] for _ in 1:nthreads]  # Intensity AC
    ac2_bootstrap_threads = [[Float64[] for _ in 1:length(lags_ac)] for _ in 1:nthreads]  # Intensity AC
    ac1ON_bootstrap_threads = [[Float64[] for _ in 1:length(lags_ac)] for _ in 1:nthreads]  # ON state AC
    ac2ON_bootstrap_threads = [[Float64[] for _ in 1:length(lags_ac)] for _ in 1:nthreads]  # ON state AC
    
    # Parallel bootstrap loop
    Threads.@threads for b in 1:n_bootstrap_cov
        thread_id = Threads.threadid()
        # Resample trials with replacement
        bootstrap_trials = StatsBase.sample(1:ntrials, ntrials, replace=true)
        
        # Compute bootstrap means from resampled trials
        bootstrap_mean1 = mean([mean1_trials[i] for i in bootstrap_trials])  # ON state mean
        bootstrap_mean2 = mean([mean2_trials[i] for i in bootstrap_trials])  # ON state mean
        bootstrap_mean1_intensity = mean([mean1_intensity_trials[i] for i in bootstrap_trials])  # Intensity mean
        bootstrap_mean2_intensity = mean([mean2_intensity_trials[i] for i in bootstrap_trials])  # Intensity mean
        
        # Compute bootstrap raw moments
        bootstrap_cc_raw = [mean([cc_raw_trials[bootstrap_trials[i]][j] for i in 1:ntrials]) for j in 1:length(lags)]  # ON state CC
        bootstrap_cc_intensity_raw = [mean([cc_intensity_raw_trials[bootstrap_trials[i]][j] for i in 1:ntrials]) for j in 1:length(lags)]  # Intensity CC
        bootstrap_ac1_raw = [mean([ac1_raw_trials[bootstrap_trials[i]][j] for i in 1:ntrials]) for j in 1:length(lags_ac)]  # Intensity AC
        bootstrap_ac2_raw = [mean([ac2_raw_trials[bootstrap_trials[i]][j] for i in 1:ntrials]) for j in 1:length(lags_ac)]  # Intensity AC
        bootstrap_ac1ON_raw = [mean([ac1ON_raw_trials[bootstrap_trials[i]][j] for i in 1:ntrials]) for j in 1:length(lags_ac)]  # ON state AC
        bootstrap_ac2ON_raw = [mean([ac2ON_raw_trials[bootstrap_trials[i]][j] for i in 1:ntrials]) for j in 1:length(lags_ac)]  # ON state AC
        
        # Compute bootstrap covariances (mean-subtracted)
        # Use GLOBAL all-trial means for centering to match main computation
        bootstrap_cc = bootstrap_cc_raw .- (mean1 * mean2)  # ON state CC (unnormalized, centered by global mean)
        # Normalize by GLOBAL all-trial means for consistency
        bootstrap_cc_normalized = bootstrap_cc ./ (mean1 * mean2)  # ON state CC (normalized by global mean)
        bootstrap_cc_intensity = bootstrap_cc_intensity_raw .- (mean1_intensity * mean2_intensity)  # Intensity CC (centered by global mean)
        bootstrap_ac1 = bootstrap_ac1_raw .- (mean1_intensity^2)  # Intensity AC1 (centered by global mean)
        bootstrap_ac2 = bootstrap_ac2_raw .- (mean2_intensity^2)  # Intensity AC2 (centered by global mean)
        bootstrap_ac1ON = bootstrap_ac1ON_raw .- (mean1^2)  # ON state AC1 (centered by global mean)
        bootstrap_ac2ON = bootstrap_ac2ON_raw .- (mean2^2)  # ON state AC2 (centered by global mean)
        
        # Store bootstrap values for each lag in thread-local arrays
        for j in 1:length(lags)
            push!(cc_bootstrap_threads[thread_id][j], bootstrap_cc_normalized[j])  # ON state CC (normalized)
            push!(cc_intensity_bootstrap_threads[thread_id][j], bootstrap_cc_intensity[j])  # Intensity CC
        end
        for j in 1:length(lags_ac)
            push!(ac1_bootstrap_threads[thread_id][j], bootstrap_ac1[j])  # Intensity AC1
            push!(ac2_bootstrap_threads[thread_id][j], bootstrap_ac2[j])  # Intensity AC2
            push!(ac1ON_bootstrap_threads[thread_id][j], bootstrap_ac1ON[j])  # ON state AC1
            push!(ac2ON_bootstrap_threads[thread_id][j], bootstrap_ac2ON[j])  # ON state AC2
        end
    end
    
    # Combine thread-local results
    cc_bootstrap = [Float64[] for _ in 1:length(lags)]  # ON state CC
    cc_intensity_bootstrap = [Float64[] for _ in 1:length(lags)]  # Intensity CC
    ac1_bootstrap = [Float64[] for _ in 1:length(lags_ac)]  # Intensity AC1
    ac2_bootstrap = [Float64[] for _ in 1:length(lags_ac)]  # Intensity AC2
    ac1ON_bootstrap = [Float64[] for _ in 1:length(lags_ac)]  # ON state AC1
    ac2ON_bootstrap = [Float64[] for _ in 1:length(lags_ac)]  # ON state AC2
    for t in 1:nthreads
        for j in 1:length(lags)
            append!(cc_bootstrap[j], cc_bootstrap_threads[t][j])  # ON state CC
            append!(cc_intensity_bootstrap[j], cc_intensity_bootstrap_threads[t][j])  # Intensity CC
        end
        for j in 1:length(lags_ac)
            append!(ac1_bootstrap[j], ac1_bootstrap_threads[t][j])  # Intensity AC1
            append!(ac2_bootstrap[j], ac2_bootstrap_threads[t][j])  # Intensity AC2
            append!(ac1ON_bootstrap[j], ac1ON_bootstrap_threads[t][j])  # ON state AC1
            append!(ac2ON_bootstrap[j], ac2ON_bootstrap_threads[t][j])  # ON state AC2
        end
    end
    
    # Compute percentiles for each lag
    # ON state CC
    cc_lower = [quantile(cc_bootstrap[j], 0.025) for j in 1:length(lags)]
    cc_median = [median(cc_bootstrap[j]) for j in 1:length(lags)]
    cc_upper = [quantile(cc_bootstrap[j], 0.975) for j in 1:length(lags)]
    # Intensity CC
    cc_intensity_lower = [quantile(cc_intensity_bootstrap[j], 0.025) for j in 1:length(lags)]
    cc_intensity_median = [median(cc_intensity_bootstrap[j]) for j in 1:length(lags)]
    cc_intensity_upper = [quantile(cc_intensity_bootstrap[j], 0.975) for j in 1:length(lags)]
    # Intensity AC
    ac1_lower = [quantile(ac1_bootstrap[j], 0.025) for j in 1:length(lags_ac)]
    ac1_median = [median(ac1_bootstrap[j]) for j in 1:length(lags_ac)]
    ac1_upper = [quantile(ac1_bootstrap[j], 0.975) for j in 1:length(lags_ac)]
    ac2_lower = [quantile(ac2_bootstrap[j], 0.025) for j in 1:length(lags_ac)]
    ac2_median = [median(ac2_bootstrap[j]) for j in 1:length(lags_ac)]
    ac2_upper = [quantile(ac2_bootstrap[j], 0.975) for j in 1:length(lags_ac)]
    # ON state AC
    ac1ON_lower = [quantile(ac1ON_bootstrap[j], 0.025) for j in 1:length(lags_ac)]
    ac1ON_median = [median(ac1ON_bootstrap[j]) for j in 1:length(lags_ac)]
    ac1ON_upper = [quantile(ac1ON_bootstrap[j], 0.975) for j in 1:length(lags_ac)]
    ac2ON_lower = [quantile(ac2ON_bootstrap[j], 0.025) for j in 1:length(lags_ac)]
    ac2ON_median = [median(ac2ON_bootstrap[j]) for j in 1:length(lags_ac)]
    ac2ON_upper = [quantile(ac2ON_bootstrap[j], 0.975) for j in 1:length(lags_ac)]
    
    # Compute norms from final mean-subtracted covariances
    # Use intensity CC for norms (cc_theory = cc, cc_mean = cc_intensity_mean)
    linf_norm = [maximum(abs.(cc .- cc_intensity_mean))]
    l2_norm = [sqrt(sum((cc .- cc_intensity_mean) .^ 2))]
    
    # Standard error of cross-covariances (variance of mean across trials)
    cc_se = sqrt.(cc_var ./ ntrials)  # ON state CC SE
    cc_intensity_se = sqrt.(cc_intensity_var ./ ntrials)  # Intensity CC SE
    
    # Normalize ON state cross-covariances (matching IDL format)
    # Normalize theoretical ccON
    ccON = ccON ./ (mON1 * mON2)  # Normalize theoretical ON state cross-covariance
    # Normalize empirical ccON
    cc_ON_mean = cc_mean ./ (mean1 * mean2)  # Normalized ON state cross-covariance (empirical)
    
    return (
        # Intensity correlations (theory: ac1, ac2, cc)
        cc_theory=cc, cc_mean=cc_intensity_mean, cc_se=cc_intensity_se, cc_lower=cc_intensity_lower, cc_median=cc_intensity_median, cc_upper=cc_intensity_upper,  # Intensity CC
        ac1_mean=ac1_mean, ac1_theory=ac1, ac1_lower=ac1_lower, ac1_median=ac1_median, ac1_upper=ac1_upper,  # Intensity AC1
        ac2_mean=ac2_mean, ac2_theory=ac2, ac2_lower=ac2_lower, ac2_median=ac2_median, ac2_upper=ac2_upper,  # Intensity AC2
        # ON state correlations (theory: ac1ON, ac2ON, ccON) - all normalized
        ccON_theory=ccON, ccON_mean=cc_ON_mean, ccON_se=cc_se ./ (mean1 * mean2), ccON_lower=cc_lower, ccON_median=cc_median, ccON_upper=cc_upper,  # ON state CC (normalized)
        cc_ON_theory=ccON, cc_ON_mean=cc_ON_mean,  # Normalized ON state CC (same as ccON, kept for compatibility)
        ac1ON_mean=ac1ON_mean, ac1ON_theory=ac1ON, ac1ON_lower=ac1ON_lower, ac1ON_median=ac1ON_median, ac1ON_upper=ac1ON_upper,  # ON state AC1
        ac2ON_mean=ac2ON_mean, ac2ON_theory=ac2ON, ac2ON_lower=ac2ON_lower, ac2ON_median=ac2ON_median, ac2ON_upper=ac2ON_upper,  # ON state AC2
        # Other statistics
        linf_norm=linf_norm, l2_norm=l2_norm, lags=lags,
        mean1=mean1, mean2=mean2, mean1_se=mean1_se, mean2_se=mean2_se, 
        mean1_lower=mean1_lower, mean1_median=mean1_median, mean1_upper=mean1_upper, 
        mean2_lower=mean2_lower, mean2_median=mean2_median, mean2_upper=mean2_upper, 
        var1=var1, var2=var2, m1=mON1, m2=mON2, v1=v1, v2=v2)
end

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
- Uses `StatsBase.crosscov` with `demean=true` which:
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
function data_covariance(traces, lags)
    result = compute_covariance(traces, lags, bootstrap=false)
    return result.ac1, result.ac2, result.cc, result.lags
end

"""
    compute_covariance(traces, lags; bootstrap=false, n_bootstrap=10000)

Core covariance computing engine. Computes cross-covariance and autocovariances from traces.

# Arguments
- `traces::Vector{Matrix{Float64}}`: Vector of trace matrices, each with columns [unit1, unit2]
- `lags::Vector{Int}`: Time lags for correlation calculation
- `bootstrap::Bool=false`: Whether to compute bootstrap confidence intervals
- `n_bootstrap::Int=10000`: Number of bootstrap replicates (if bootstrap=true)

# Returns
Named tuple with:
- `cc`: Cross-covariance (E[x(t)y(t+τ)] - E[x]E[y])
- `ac1`: Autocovariance for unit1 (E[x(t)x(t+τ)] - E[x]^2)
- `ac2`: Autocovariance for unit2 (E[y(t)y(t+τ)] - E[y]^2)
- `mean1`, `mean2`: Mean values
- `lags`: Time lags

If `bootstrap=true`, also includes confidence intervals for all quantities.

# Algorithm
Uses `StatsBase.crosscov` with `demean=true` to compute demeaned covariances.
Averages across all traces.
"""
function compute_covariance(traces::Vector{Matrix{Float64}}, lags::Vector{Int}; bootstrap=false, n_bootstrap=10000)
    n_traces = length(traces)
    
    if n_traces == 0
        error("No traces provided")
    end
    
    # Compute overall means from all traces
    trace1_all = Float64[]
    trace2_all = Float64[]
    for t in traces
        append!(trace1_all, t[:, 1])
        append!(trace2_all, t[:, 2])
    end
    mean1 = mean(trace1_all)
    mean2 = mean(trace2_all)
    
    # Compute cross-covariance and autocovariances using StatsBase.crosscov
    cc_sum = zeros(length(lags))
    ac1_sum = zeros(length(lags))
    ac2_sum = zeros(length(lags))
    
    for t in traces
        # Use StatsBase.crosscov with demean=true
        cc_sum += StatsBase.crosscov(t[:, 1], t[:, 2], lags, demean=true)
        ac1_sum += StatsBase.crosscov(t[:, 1], t[:, 1], lags, demean=true)
        ac2_sum += StatsBase.crosscov(t[:, 2], t[:, 2], lags, demean=true)
    end
    
    # Average across traces
    cc = cc_sum / n_traces
    ac1 = ac1_sum / n_traces
    ac2 = ac2_sum / n_traces
    
    result = (
        cc=cc,
        ac1=ac1,
        ac2=ac2,
        mean1=mean1,
        mean2=mean2,
        lags=lags
    )
    
    # Bootstrap if requested
    if bootstrap
        cc_bootstrap = [Float64[] for _ in 1:length(lags)]
        ac1_bootstrap = [Float64[] for _ in 1:length(lags)]
        ac2_bootstrap = [Float64[] for _ in 1:length(lags)]
        mean1_bootstrap = Float64[]
        mean2_bootstrap = Float64[]
        
        for b in 1:n_bootstrap
            # Resample traces with replacement
            bootstrap_traces = StatsBase.sample(traces, n_traces, replace=true)
            
            # Compute bootstrap means
            bootstrap_trace1_all = Float64[]
            bootstrap_trace2_all = Float64[]
            for t in bootstrap_traces
                append!(bootstrap_trace1_all, t[:, 1])
                append!(bootstrap_trace2_all, t[:, 2])
            end
            push!(mean1_bootstrap, mean(bootstrap_trace1_all))
            push!(mean2_bootstrap, mean(bootstrap_trace2_all))
            
            # Compute bootstrap correlations
            bootstrap_cc_sum = zeros(length(lags))
            bootstrap_ac1_sum = zeros(length(lags))
            bootstrap_ac2_sum = zeros(length(lags))
            
            for t in bootstrap_traces
                bootstrap_cc_sum += StatsBase.crosscov(t[:, 1], t[:, 2], lags, demean=true)
                bootstrap_ac1_sum += StatsBase.crosscov(t[:, 1], t[:, 1], lags, demean=true)
                bootstrap_ac2_sum += StatsBase.crosscov(t[:, 2], t[:, 2], lags, demean=true)
            end
            
            bootstrap_cc = bootstrap_cc_sum / n_traces
            bootstrap_ac1 = bootstrap_ac1_sum / n_traces
            bootstrap_ac2 = bootstrap_ac2_sum / n_traces
            
            for j in 1:length(lags)
                push!(cc_bootstrap[j], bootstrap_cc[j])
                push!(ac1_bootstrap[j], bootstrap_ac1[j])
                push!(ac2_bootstrap[j], bootstrap_ac2[j])
            end
        end
        
        # Compute percentiles
        cc_lower = [quantile(cc_bootstrap[j], 0.025) for j in 1:length(lags)]
        cc_median = [median(cc_bootstrap[j]) for j in 1:length(lags)]
        cc_upper = [quantile(cc_bootstrap[j], 0.975) for j in 1:length(lags)]
        
        ac1_lower = [quantile(ac1_bootstrap[j], 0.025) for j in 1:length(lags)]
        ac1_median = [median(ac1_bootstrap[j]) for j in 1:length(lags)]
        ac1_upper = [quantile(ac1_bootstrap[j], 0.975) for j in 1:length(lags)]
        
        ac2_lower = [quantile(ac2_bootstrap[j], 0.025) for j in 1:length(lags)]
        ac2_median = [median(ac2_bootstrap[j]) for j in 1:length(lags)]
        ac2_upper = [quantile(ac2_bootstrap[j], 0.975) for j in 1:length(lags)]
        
        mean1_lower = quantile(mean1_bootstrap, 0.025)
        mean1_median = median(mean1_bootstrap)
        mean1_upper = quantile(mean1_bootstrap, 0.975)
        
        mean2_lower = quantile(mean2_bootstrap, 0.025)
        mean2_median = median(mean2_bootstrap)
        mean2_upper = quantile(mean2_bootstrap, 0.975)
        
        return merge(result, (
            cc_lower=cc_lower,
            cc_median=cc_median,
            cc_upper=cc_upper,
            ac1_lower=ac1_lower,
            ac1_median=ac1_median,
            ac1_upper=ac1_upper,
            ac2_lower=ac2_lower,
            ac2_median=ac2_median,
            ac2_upper=ac2_upper,
            mean1_lower=mean1_lower,
            mean1_median=mean1_median,
            mean1_upper=mean1_upper,
            mean2_lower=mean2_lower,
            mean2_median=mean2_median,
            mean2_upper=mean2_upper
        ))
    end
    
    return result
end

# Backward compatibility wrapper
function data_covariance(traces, lags)
    result = compute_covariance(traces, lags, bootstrap=false)
    return result.ac1, result.ac2, result.cc, result.lags
end

"""
    crosscovariance_gof_test(data_traces, r, transitions, G, R, S, insertstep, coupling, interval, probfn, lags; ntrials=1000, trial_time=720.0)

Simulation-based goodness-of-fit test comparing empirical cross-covariance to model predictions.

This function performs a comprehensive test by:
1. Computing empirical cross-covariance from data traces
2. Computing predicted cross-covariance from model parameters
3. Simulating many traces from the model and computing their cross-covariances
4. Comparing where the empirical cross-covariance falls in the simulation distribution

# Arguments
- `data_traces::Vector`: Vector of empirical trace data (matrices with columns [enhancer, gene])
- `r::Vector`: Fitted rate parameters from the model
- `transitions::Tuple`: Model transition structure
- `G, R, S, insertstep`: Model structure parameters
- `coupling::Tuple`: Coupling structure for coupled models
- `interval::Float64`: Time interval between trace points
- `probfn::Function`: Probability function for noise (e.g., `prob_Gaussian`)
- `lags::Vector{Int}`: Time lags for cross-covariance calculation
- `ntrials::Int=1000`: Number of simulation trials (default: 1000)
- `trial_time::Float64=720.0`: Length of each simulated trace (default: 720.0 minutes)

# Returns
- `NamedTuple` with fields:
  - `empirical_cc::Vector`: Empirical cross-covariance from data
  - `predicted_cc::Vector`: Predicted cross-covariance from model
  - `simulation_cc_mean::Vector`: Mean cross-covariance across simulations
  - `simulation_cc_std::Vector`: Standard deviation of cross-covariance across simulations
  - `simulation_cc_samples::Matrix`: All simulation cross-covariances (ntrials × length(lags))
  - `lags::Vector`: Time lags
  - `l2_norm::Float64`: L² norm of difference between empirical and predicted
  - `linf_norm::Float64`: L∞ norm (max absolute difference)
  - `pearson_corr::Float64`: Pearson correlation between empirical and predicted
  - `percentile_rank::Float64`: Percentile rank of empirical L² norm in simulation distribution
  - `p_value::Float64`: Two-tailed p-value (proportion of simulations with larger L² norm)

# Notes
- Uses `data_covariance` to compute empirical cross-covariance (handles edge effects and trends)
- Uses `covariance_functions` to compute predicted cross-covariance
- Simulates traces using `simulate_trace_vector` and computes cross-covariances
- Accounts for finite-sample effects by comparing to simulation distribution
- Edge effects and trends are automatically handled by `StatsBase.crosscov` with `demean=true`

# Examples
```julia
# Perform goodness-of-fit test
lags = collect(-60:1:60)
results = crosscovariance_gof_test(
    data_traces, fitted_rates, transitions, G, R, S, insertstep,
    coupling, interval, prob_Gaussian, lags;
    ntrials=1000, trial_time=720.0
)

# Check p-value
println("P-value: ", results.p_value)
println("Percentile rank: ", results.percentile_rank)

# Visualize: empirical vs predicted vs simulation distribution
using Plots
using Printf
plot(results.lags, results.empirical_cc, label="Empirical")
plot!(results.lags, results.predicted_cc, label="Predicted")
plot!(results.lags, results.simulation_cc_mean, ribbon=results.simulation_cc_std, 
      label="Simulation mean ± std", alpha=0.3)
```
"""
function crosscovariance_gof_test(data_traces, r, transitions, G, R, S, insertstep, coupling, interval, probfn, lags; ntrials=1000, trial_time=720.0)
    # 1. Compute empirical cross-covariance from data
    _, _, empirical_cc, _ = data_covariance(data_traces, lags)
    
    # 2. Compute predicted cross-covariance from model
    _, _, predicted_cc, _, _, _, _, _, _, _ = covariance_functions(r, transitions, G, R, S, insertstep, interval, probfn, coupling, lags[lags.>=0])
    # Extend to negative lags (symmetric)
    if any(lags.<0)
        lags_positive = lags[lags.>=0]
        predicted_cc_full = vcat(reverse(predicted_cc[2:end]), predicted_cc)
    else
        predicted_cc_full = predicted_cc
    end
    
    # 3. Simulate many traces and compute their cross-covariances
    simulation_cc_samples = zeros(ntrials, length(lags))
    l2_norms = zeros(ntrials)
    
    for i in 1:ntrials
        # Simulate a trace from the model
        sim_trace = simulate_trace_vector(r, transitions, G, R, S, insertstep, coupling, interval, trial_time, 1, col=2)
        
        # Compute cross-covariance from this simulation
        _, _, sim_cc, _ = data_covariance(sim_trace, lags)
        simulation_cc_samples[i, :] = sim_cc
        
        # Compute L² norm for this simulation
        l2_norms[i] = sqrt(sum((sim_cc .- predicted_cc_full).^2))
    end
    
    # 4. Compute statistics
    simulation_cc_mean = vec(Statistics.mean(simulation_cc_samples, dims=1))
    simulation_cc_std = vec(Statistics.std(simulation_cc_samples, dims=1))
    
    # Compare empirical to predicted
    l2_norm = sqrt(sum((empirical_cc .- predicted_cc_full).^2))
    linf_norm = maximum(abs.(empirical_cc .- predicted_cc_full))
    
    # Pearson correlation
    pearson_corr = StatsBase.cor(empirical_cc, predicted_cc_full)
    
    # Percentile rank: where does empirical L² norm fall in simulation distribution?
    percentile_rank = sum(l2_norms .<= l2_norm) / ntrials
    
    # Two-tailed p-value: proportion of simulations with larger or equal L² norm
    # (more extreme than observed)
    p_value = 2 * min(percentile_rank, 1.0 - percentile_rank)
    
    return (
        empirical_cc=empirical_cc,
        predicted_cc=predicted_cc_full,
        simulation_cc_mean=simulation_cc_mean,
        simulation_cc_std=simulation_cc_std,
        simulation_cc_samples=simulation_cc_samples,
        lags=lags,
        l2_norm=l2_norm,
        linf_norm=linf_norm,
        pearson_corr=pearson_corr,
        percentile_rank=percentile_rank,
        p_value=p_value
    )
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
    getratefile(files, gene)

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
- Uses `StatsBase.crosscov` with `demean=true` which:
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

"""
    bootstrap_crosscovariance(data_traces, lags; n_bootstrap=10000)

Compute bootstrap confidence intervals for cross-covariance from empirical trace data.

# Arguments
- `data_traces::Vector{Matrix}`: Vector of trace matrices, each with columns [enhancer, gene]
  - 137 traces for simultaneous measurements
- `lags::Vector{Int}`: Time lags for cross-covariance calculation
- `n_bootstrap::Int=10000`: Number of bootstrap replicates

# Returns
Named tuple with:
- `empirical_cc`: Mean empirical cross-covariance across all traces
- `empirical_cc_lower`: Lower bootstrap CI (2.5th percentile)
- `empirical_cc_median`: Median bootstrap value
- `empirical_cc_upper`: Upper bootstrap CI (97.5th percentile)
- `empirical_cc_std`: Standard deviation across bootstrap samples
- `lags`: Time lags

# Algorithm
Matches `simulate_trials` approach:
1. Compute raw cross-correlation with `demean=false`
2. Compute overall means from all traces
3. Subtract `mean1 * mean2` to get cross-covariance
4. Bootstrap resample traces with replacement
5. Compute cross-covariance for each bootstrap sample
6. Compute percentiles for confidence intervals

# Notes
- Maintains trace pairing (enhancer and gene from same trace stay together)
- Uses thread-local storage for parallel bootstrap computation
"""
function bootstrap_crosscovariance(data_traces::Vector{Matrix{Float64}}, lags; n_bootstrap=10000)
    n_traces = length(data_traces)
    
    # Compute overall means from all traces (for consistency with simulate_trials)
    trace1_all = Float64[]
    trace2_all = Float64[]
    for t in data_traces
        append!(trace1_all, t[:, 1])
        append!(trace2_all, t[:, 2])
    end
    mean1 = mean(trace1_all)
    mean2 = mean(trace2_all)
    
    # Compute empirical cross-covariance from all data
    cc_raw_all = zeros(length(lags))
    for t in data_traces
        cc_raw_all .+= StatsBase.crosscov(t[:, 1], t[:, 2], lags, demean=false)
    end
    cc_raw_mean = cc_raw_all / n_traces
    empirical_cc = cc_raw_mean .- (mean1 * mean2)
    
    # Bootstrap resampling
    nthreads = Threads.nthreads()
    cc_bootstrap_threads = [[Float64[] for _ in 1:length(lags)] for _ in 1:nthreads]
    
    Threads.@threads for b in 1:n_bootstrap
        thread_id = Threads.threadid()
        # Resample traces with replacement (maintains pairing)
        bootstrap_traces = StatsBase.sample(1:n_traces, n_traces, replace=true)
        
        # Compute raw cross-correlation for bootstrap sample
        cc_raw_bootstrap = zeros(length(lags))
        mean1_bootstrap = 0.0
        mean2_bootstrap = 0.0
        n_points_bootstrap = 0
        
        for idx in bootstrap_traces
            t = data_traces[idx]
            cc_raw_bootstrap .+= StatsBase.crosscov(t[:, 1], t[:, 2], lags, demean=false)
            mean1_bootstrap += sum(t[:, 1])
            mean2_bootstrap += sum(t[:, 2])
            n_points_bootstrap += size(t, 1)
        end
        
        # Compute means for bootstrap sample
        mean1_bs = mean1_bootstrap / n_points_bootstrap
        mean2_bs = mean2_bootstrap / n_points_bootstrap
        
        # Compute cross-covariance (subtract means)
        cc_bootstrap = (cc_raw_bootstrap / n_traces) .- (mean1_bs * mean2_bs)
        
        # Store in thread-local arrays
        for j in 1:length(lags)
            push!(cc_bootstrap_threads[thread_id][j], cc_bootstrap[j])
        end
    end
    
    # Combine thread-local arrays
    cc_bootstrap = [Float64[] for _ in 1:length(lags)]
    for t in 1:nthreads
        for j in 1:length(lags)
            append!(cc_bootstrap[j], cc_bootstrap_threads[t][j])
        end
    end
    
    # Compute percentiles
    empirical_cc_lower = [quantile(cc_bootstrap[j], 0.025) for j in 1:length(lags)]
    empirical_cc_median = [median(cc_bootstrap[j]) for j in 1:length(lags)]
    empirical_cc_upper = [quantile(cc_bootstrap[j], 0.975) for j in 1:length(lags)]
    empirical_cc_std = [std(cc_bootstrap[j]) for j in 1:length(lags)]
    
    return (
        empirical_cc=empirical_cc,
        empirical_cc_lower=empirical_cc_lower,
        empirical_cc_median=empirical_cc_median,
        empirical_cc_upper=empirical_cc_upper,
        empirical_cc_std=empirical_cc_std,
        lags=lags
    )
end

"""
    score_crosscovariance_model(data_traces, r, transitions, G, R, S, insertstep, coupling, interval, probfn, lags; n_bootstrap=10000)

Score a fitted model's cross-covariance predictions against empirical data using bootstrap confidence intervals.

# Arguments
- `data_traces::Vector{Matrix}`: Vector of 137 trace matrices [enhancer, gene]
- `r::Vector{Float64}`: Fitted rate parameters
- `transitions`, `G`, `R`, `S`, `insertstep`, `coupling`: Model structure
- `interval::Float64`: Time interval between frames
- `probfn`: Observation distribution function
- `lags::Vector{Int}`: Time lags for comparison
- `n_bootstrap::Int=10000`: Number of bootstrap replicates

# Returns
Named tuple with scoring metrics:
- `theoretical_cc`: Theoretical cross-covariance from model
- `empirical_cc`: Mean empirical cross-covariance
- `empirical_cc_lower`, `empirical_cc_upper`: Bootstrap confidence intervals
- `l2_norm`: L² distance between theory and empirical
- `linf_norm`: L∞ distance (maximum absolute difference)
- `z_scores`: Per-lag z-scores: (theory - empirical_mean) / empirical_std
- `coverage`: Boolean vector indicating if theory is within CI for each lag
- `coverage_fraction`: Fraction of lags where theory is within CI
- `lags`: Time lags
- `m1`, `m2`: Theoretical means

# Algorithm Comparison:
- **IDL**: Computes normalized cross-covariance (divides by product of means), stores unnormalized, then normalizes in bootstrap
- **Julia**: Computes unnormalized cross-covariance directly (matches theory)
- Both should give similar results after accounting for normalization differences
"""
function score_crosscovariance_model(data_traces::Vector{Matrix{Float64}}, r, transitions, G, R, S, insertstep, coupling, interval, probfn, lags; n_bootstrap=10000)
    # 1. Compute theoretical cross-covariance
    positive_lags = lags[lags .>= 0]
    if isempty(positive_lags)
        error("lags must contain at least one non-negative lag")
    end
    
    ac1_theory, ac2_theory, cc_theory, _, full_lags, m1, m2, v1, v2 = 
        StochasticGene.covariance_functions(r, transitions, G, R, S, insertstep, interval, probfn, coupling, positive_lags)
    
    # Ensure lags match - handle negative lags by using symmetry
    if any(lags .< 0)
        # For cross-covariance, we need both positive and negative lags
        # Theory gives us positive lags, we need to construct full symmetric version
        all_positive_lags = sort(unique([abs.(lags); positive_lags]))
        ac1_theory_full, ac2_theory_full, cc_theory_full, _, full_lags_full, m1, m2, v1, v2 = 
            StochasticGene.covariance_functions(r, transitions, G, R, S, insertstep, interval, probfn, coupling, all_positive_lags)
        
        # Construct symmetric cross-covariance (cc(-τ) = cc(τ) for autocovariance, but cross-covariance is not symmetric)
        # For cross-covariance, we need cc12(τ) and cc21(τ) = cc12(-τ)
        # The theory function returns cc12 for positive lags
        # We need to construct the full symmetric array
        cc_theory_symmetric = zeros(length(lags))
        for (i, lag) in enumerate(lags)
            if lag >= 0
                idx = findfirst(==(lag), full_lags_full)
                if !isnothing(idx)
                    cc_theory_symmetric[i] = cc_theory_full[idx]
                end
            else
                # For negative lags, find the corresponding positive lag
                idx = findfirst(==(abs(lag)), full_lags_full)
                if !isnothing(idx)
                    # Cross-covariance: cc12(-τ) = cc21(τ), but we computed cc12
                    # Need to check if theory returns cc12 or combined cc
                    # Based on covariance_functions, it returns combined cc with negative lags first
                    # So we need to find the corresponding entry
                    cc_theory_symmetric[i] = cc_theory_full[idx]
                end
            end
        end
        cc_theory = cc_theory_symmetric
        full_lags = full_lags_full
    end
    
    # Match lags if they don't align exactly
    if full_lags != lags
        cc_theory_matched = zeros(length(lags))
        for (i, lag) in enumerate(lags)
            idx = findfirst(==(lag), full_lags)
            if !isnothing(idx)
                cc_theory_matched[i] = cc_theory[idx]
            else
                # Linear interpolation if lag not found
                idx1 = findlast(<(lag), full_lags)
                idx2 = findfirst(>(lag), full_lags)
                if !isnothing(idx1) && !isnothing(idx2)
                    w1 = (full_lags[idx2] - lag) / (full_lags[idx2] - full_lags[idx1])
                    w2 = (lag - full_lags[idx1]) / (full_lags[idx2] - full_lags[idx1])
                    cc_theory_matched[i] = w1 * cc_theory[idx1] + w2 * cc_theory[idx2]
                end
            end
        end
        cc_theory = cc_theory_matched
    end
    
    # 2. Compute bootstrap empirical cross-covariance
    bootstrap_result = bootstrap_crosscovariance(data_traces, lags; n_bootstrap=n_bootstrap)
    
    # 3. Compute scoring metrics
    l2_norm = sqrt(sum((cc_theory .- bootstrap_result.empirical_cc).^2))
    linf_norm = maximum(abs.(cc_theory .- bootstrap_result.empirical_cc))
    
    # Compute z-scores: (theory - empirical_mean) / empirical_std
    z_scores = (cc_theory .- bootstrap_result.empirical_cc) ./ bootstrap_result.empirical_cc_std
    z_scores[isnan.(z_scores)] .= 0.0  # Handle division by zero
    
    # Coverage: whether theory is within bootstrap CI
    coverage = (cc_theory .>= bootstrap_result.empirical_cc_lower) .& 
               (cc_theory .<= bootstrap_result.empirical_cc_upper)
    
    return (
        theoretical_cc=cc_theory,
        empirical_cc=bootstrap_result.empirical_cc,
        empirical_cc_lower=bootstrap_result.empirical_cc_lower,
        empirical_cc_median=bootstrap_result.empirical_cc_median,
        empirical_cc_upper=bootstrap_result.empirical_cc_upper,
        empirical_cc_std=bootstrap_result.empirical_cc_std,
        l2_norm=l2_norm,
        linf_norm=linf_norm,
        z_scores=z_scores,
        coverage=coverage,
        coverage_fraction=mean(coverage),
        lags=lags,
        m1=m1,
        m2=m2
    )
end

"""
    compute_normalized_crosscovariance(fitted_traces, lags; bootstrap=false, n_bootstrap=10000)

Compute normalized cross-covariance (CC_ON format) from precomputed fitted traces.

# Arguments
- `fitted_traces::Vector{Matrix{Float64}}`: Vector of fitted trace matrices, each with columns [enhancer, gene]
- `lags::Vector{Int}`: Time lags for cross-covariance calculation
- `bootstrap::Bool=false`: Whether to compute bootstrap confidence intervals
- `n_bootstrap::Int=10000`: Number of bootstrap replicates (if bootstrap=true)

# Returns
If `bootstrap=false`: Named tuple with:
- `cc_normalized`: Normalized cross-covariance: (mean(xy) - mean(x)*mean(y)) / (mean(x)*mean(y))
- `cc_unnormalized`: Unnormalized cross-covariance: mean(xy) - mean(x)*mean(y)
- `mean1`, `mean2`: Means of enhancer and gene traces
- `lags`: Time lags

If `bootstrap=true`: Also includes:
- `cc_normalized_lower`, `cc_normalized_upper`: Bootstrap CI for normalized CC
- `cc_normalized_median`: Median bootstrap value
- `cc_normalized_std`: Standard deviation across bootstrap samples

# Algorithm
Matches IDL algorithm (line 709):
- Computes: `Gn = (n * Go) / (Mdirect * Mdelayed) - 1.0`
- Where: `Go = sum(xy)`, `Mdirect = sum(x)`, `Mdelayed = sum(y)`, `n = number of pairs`
- This gives: `Gn = mean(xy) / (mean(x) * mean(y)) - 1.0 = (mean(xy) - mean(x)*mean(y)) / (mean(x)*mean(y))`

This is the normalized cross-covariance format used by experimentalists (CC_ON).

# Notes
- Normalized format allows comparison across different mean intensity levels
- Bootstrap resampling maintains trace pairing (enhancer and gene from same trace stay together)
"""
function compute_normalized_crosscovariance(fitted_traces::Vector{Matrix{Float64}}, lags; bootstrap=false, n_bootstrap=10000)
    n_traces = length(fitted_traces)
    
    # Compute overall means from all traces
    trace1_all = Float64[]
    trace2_all = Float64[]
    for t in fitted_traces
        append!(trace1_all, t[:, 1])
        append!(trace2_all, t[:, 2])
    end
    mean1 = mean(trace1_all)
    mean2 = mean(trace2_all)
    
    # Compute raw cross-correlation (unnormalized)
    cc_raw_all = zeros(length(lags))
    for t in fitted_traces
        cc_raw_all .+= StatsBase.crosscov(t[:, 1], t[:, 2], lags, demean=false)
    end
    cc_raw_mean = cc_raw_all / n_traces
    
    # Compute unnormalized cross-covariance
    cc_unnormalized = cc_raw_mean .- (mean1 * mean2)
    
    # Compute normalized cross-covariance (CC_ON format)
    # Normalize by product of means: (mean(xy) - mean(x)*mean(y)) / (mean(x)*mean(y))
    cc_normalized = cc_unnormalized ./ (mean1 * mean2)
    
    result = (
        cc_normalized=cc_normalized,
        cc_unnormalized=cc_unnormalized,
        mean1=mean1,
        mean2=mean2,
        lags=lags
    )
    
    if !bootstrap
        return result
    end
    
    # Bootstrap resampling for confidence intervals
    nthreads = Threads.nthreads()
    cc_normalized_bootstrap_threads = [[Float64[] for _ in 1:length(lags)] for _ in 1:nthreads]
    
    Threads.@threads for b in 1:n_bootstrap
        thread_id = Threads.threadid()
        # Resample traces with replacement (maintains pairing)
        bootstrap_traces = StatsBase.sample(1:n_traces, n_traces, replace=true)
        
        # Compute means for bootstrap sample
        mean1_bs = 0.0
        mean2_bs = 0.0
        n_points_bs = 0
        cc_raw_bs = zeros(length(lags))
        
        for idx in bootstrap_traces
            t = fitted_traces[idx]
            cc_raw_bs .+= StatsBase.crosscov(t[:, 1], t[:, 2], lags, demean=false)
            mean1_bs += sum(t[:, 1])
            mean2_bs += sum(t[:, 2])
            n_points_bs += size(t, 1)
        end
        
        mean1_bs /= n_points_bs
        mean2_bs /= n_points_bs
        cc_raw_mean_bs = cc_raw_bs / n_traces
        
        # Compute normalized cross-covariance for bootstrap sample
        cc_unnormalized_bs = cc_raw_mean_bs .- (mean1_bs * mean2_bs)
        cc_normalized_bs = cc_unnormalized_bs ./ (mean1_bs * mean2_bs)
        
        # Store in thread-local arrays
        for j in 1:length(lags)
            push!(cc_normalized_bootstrap_threads[thread_id][j], cc_normalized_bs[j])
        end
    end
    
    # Combine thread-local arrays
    cc_normalized_bootstrap = [Float64[] for _ in 1:length(lags)]
    for t in 1:nthreads
        for j in 1:length(lags)
            append!(cc_normalized_bootstrap[j], cc_normalized_bootstrap_threads[t][j])
        end
    end
    
    # Compute percentiles
    cc_normalized_lower = [quantile(cc_normalized_bootstrap[j], 0.025) for j in 1:length(lags)]
    cc_normalized_median = [median(cc_normalized_bootstrap[j]) for j in 1:length(lags)]
    cc_normalized_upper = [quantile(cc_normalized_bootstrap[j], 0.975) for j in 1:length(lags)]
    cc_normalized_std = [std(cc_normalized_bootstrap[j]) for j in 1:length(lags)]
    
    return (
        cc_normalized=cc_normalized,
        cc_unnormalized=cc_unnormalized,
        cc_normalized_lower=cc_normalized_lower,
        cc_normalized_median=cc_normalized_median,
        cc_normalized_upper=cc_normalized_upper,
        cc_normalized_std=cc_normalized_std,
        mean1=mean1,
        mean2=mean2,
        lags=lags
    )
end

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
    load_raw_intensity_csv(unit1_filename, unit2_filename; intensity_col=nothing, trace_id_col=nothing)

Load raw intensity traces from CSV files.

# Arguments
- `unit1_filename::String`: Path to CSV file for unit 1 (enhancer)
- `unit2_filename::String`: Path to CSV file for unit 2 (gene)
- `intensity_col::Union{String, Nothing}=nothing`: Column name for intensity (default: auto-detect "model_mean" or "Intensity" columns)
- `trace_id_col::Union{String, Nothing}=nothing`: Column name for trace ID (if multiple traces)

# Returns
- `Vector{Matrix{Float64}}`: Vector of raw intensity matrices, each with columns [unit1_intensity, unit2_intensity]
"""
function load_raw_intensity_csv(unit1_filename, unit2_filename; intensity_col=nothing, trace_id_col=nothing)
    df1 = CSV.read(unit1_filename, DataFrame)
    df2 = CSV.read(unit2_filename, DataFrame)
    
    # Auto-detect intensity columns
    if isnothing(intensity_col)
        # Look for "model_mean1", "model_mean2", etc. or "Intensity1", "Intensity2", etc.
        intensity_cols1 = filter(n -> occursin("mean", lowercase(string(n))) || occursin("intensity", lowercase(string(n))), names(df1))
        intensity_cols2 = filter(n -> occursin("mean", lowercase(string(n))) || occursin("intensity", lowercase(string(n))), names(df2))
        
        if isempty(intensity_cols1) || isempty(intensity_cols2)
            error("Could not find intensity columns in files. Found columns: unit1=$(names(df1)), unit2=$(names(df2))")
        end
        
        intensity_cols1 = sort(intensity_cols1)
        intensity_cols2 = sort(intensity_cols2)
        
        if length(intensity_cols1) != length(intensity_cols2)
            error("Different number of intensity columns: unit1 has $(length(intensity_cols1)), unit2 has $(length(intensity_cols2))")
        end
        
        traces = Matrix{Float64}[]
        for i in 1:length(intensity_cols1)
            unit1_intensity = df1[:, intensity_cols1[i]]
            unit2_intensity = df2[:, intensity_cols2[i]]
            
            # Handle missing values: replace with 0
            unit1_intensity = coalesce.(unit1_intensity, 0.0)
            unit2_intensity = coalesce.(unit2_intensity, 0.0)
            
            if length(unit1_intensity) != length(unit2_intensity)
                error("Trace $i has different lengths: unit1=$(length(unit1_intensity)), unit2=$(length(unit2_intensity))")
            end
            
            trace_matrix = hcat(unit1_intensity, unit2_intensity)
            push!(traces, trace_matrix)
        end
        
        return traces
    end
    
    # Single intensity column specified
    unit1_intensity = df1[:, intensity_col]
    unit2_intensity = df2[:, intensity_col]
    
    # Handle missing values
    unit1_intensity = coalesce.(unit1_intensity, 0.0)
    unit2_intensity = coalesce.(unit2_intensity, 0.0)
    
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
                trace_matrix = hcat(unit1_intensity[idx1], unit2_intensity[idx2])
                push!(traces, trace_matrix)
            end
            
            return traces
        else
            error("trace_id_col '$trace_id_col' not found in one or both dataframes")
        end
    else
        # Single trace
        if length(unit1_intensity) != length(unit2_intensity)
            error("Unit1 and unit2 files have different numbers of rows")
        end
        trace_matrix = hcat(unit1_intensity, unit2_intensity)
        return [trace_matrix]
    end
end

"""
    compute_binary_correlations(binary_traces, lags; bootstrap=false, n_bootstrap=10000)

Compute normalized cross-covariance and autocovariance on binary ON/OFF traces.

# Arguments
- `binary_traces::Vector{Matrix{Float64}}`: Vector of binary trace matrices [enhancer_ONOFF, gene_ONOFF]
  - Values should be 0 (OFF) or 1 (ON)
- `lags::Vector{Int}`: Time lags for correlation calculation
- `bootstrap::Bool=false`: Whether to compute bootstrap confidence intervals
- `n_bootstrap::Int=10000`: Number of bootstrap replicates (if bootstrap=true)

# Returns
Named tuple with:
- `cc_normalized`: Normalized cross-covariance (CC_ON format)
- `ac1_normalized`: Normalized autocovariance for enhancer
- `ac2_normalized`: Normalized autocovariance for gene
- `cc_unnormalized`: Unnormalized cross-covariance
- `ac1_unnormalized`: Unnormalized autocovariance for enhancer
- `ac2_unnormalized`: Unnormalized autocovariance for gene
- `mean1`, `mean2`: Mean ON probabilities (fraction ON)
- `lags`: Time lags

If `bootstrap=true`, also includes confidence intervals for all quantities.

# Algorithm
For binary traces, the normalized correlations are:
- Cross-covariance: `(P(enhancer=ON, gene=ON) - P(enhancer=ON)*P(gene=ON)) / (P(enhancer=ON)*P(gene=ON))`
- Autocovariance: `(P(state(t)=ON, state(t+τ)=ON) - P(ON)^2) / P(ON)^2`

This matches the CC_ON format used by experimentalists.
"""
function compute_binary_correlations(binary_traces::Vector{Matrix{Float64}}, lags; bootstrap=false, n_bootstrap=10000)
    n_traces = length(binary_traces)
    
    # Compute overall means (fraction ON) from all traces
    trace1_all = Float64[]
    trace2_all = Float64[]
    for t in binary_traces
        append!(trace1_all, t[:, 1])
        append!(trace2_all, t[:, 2])
    end
    mean1 = mean(trace1_all)  # P(enhancer=ON)
    mean2 = mean(trace2_all)  # P(gene=ON)
    
    # Compute raw cross-correlation and autocorrelations manually
    # Track sums for each lag to compute E[x(t)y(t+lag)] etc.
    # This matches the IDL approach: compute per trace, then average
    cc_sum_xy = zeros(length(lags))  # sum of x*y for each lag
    cc_n_pairs = zeros(Int, length(lags))  # number of pairs for each lag
    
    ac1_sum_xx = zeros(length(lags))  # sum of x*x for each lag
    ac1_n_pairs = zeros(Int, length(lags))
    
    ac2_sum_yy = zeros(length(lags))  # sum of y*y for each lag
    ac2_n_pairs = zeros(Int, length(lags))
    
    for t in binary_traces
        trace_len = size(t, 1)
        x = t[:, 1]
        y = t[:, 2]
        
        # Compute correlations manually for each lag
        # Use symmetric time range: for lag ±τ, use time points where both x(t) and y(t±τ) are valid
        # This ensures symmetry: same number of pairs for +τ and -τ
        for (i, lag) in enumerate(lags)
            # Check if this lag is valid for this trace
            lag_abs = abs(lag)
            if lag_abs >= trace_len
                continue
            end
            
            # Use symmetric range: start at max(1, 1+lag_abs) and end at min(trace_len, trace_len-lag_abs)
            # For positive lag: use t from 1 to trace_len-lag
            # For negative lag: use t from lag_abs+1 to trace_len
            # To ensure symmetry, use the intersection: t from lag_abs+1 to trace_len-lag_abs
            start_idx = lag_abs + 1
            end_idx = trace_len - lag_abs
            
            if start_idx > end_idx
                continue
            end
            
            n_pairs = end_idx - start_idx + 1
            
            if lag >= 0
                # Positive lag: E[x(t)y(t+lag)]
                # x_seg = x[t] for t from start_idx to end_idx
                # y_seg = y[t+lag] for t from start_idx to end_idx
                x_seg = x[start_idx:end_idx]
                y_seg = y[(start_idx+lag):(end_idx+lag)]
            else
                # Negative lag: E[x(t)y(t-|lag|)]
                # x_seg = x[t] for t from start_idx to end_idx
                # y_seg = y[t-lag_abs] for t from start_idx to end_idx
                x_seg = x[start_idx:end_idx]
                y_seg = y[(start_idx-lag_abs):(end_idx-lag_abs)]
            end
            
            # Verify lengths match
            if length(x_seg) != length(y_seg) || length(x_seg) != n_pairs
                error("Length mismatch: x_seg=$(length(x_seg)), y_seg=$(length(y_seg)), expected=$n_pairs, lag=$lag, trace_len=$trace_len")
            end
            
            cc_sum_xy[i] += sum(x_seg .* y_seg)
            cc_n_pairs[i] += n_pairs
            
            ac1_sum_xx[i] += sum(x_seg .* x_seg)
            ac1_n_pairs[i] += n_pairs
            
            ac2_sum_yy[i] += sum(y_seg .* y_seg)
            ac2_n_pairs[i] += n_pairs
        end
    end
    
    # Compute mean correlations (E[x(t)y(t+lag)] etc.) for each lag
    # This gives mean(xy) averaged across all traces and pairs
    cc_mean_xy = [n > 0 ? sum_xy / n : 0.0 for (n, sum_xy) in zip(cc_n_pairs, cc_sum_xy)]
    ac1_mean_xx = [n > 0 ? sum_xx / n : 0.0 for (n, sum_xx) in zip(ac1_n_pairs, ac1_sum_xx)]
    ac2_mean_yy = [n > 0 ? sum_yy / n : 0.0 for (n, sum_yy) in zip(ac2_n_pairs, ac2_sum_yy)]
    
    # Compute unnormalized covariances (using global means)
    # This matches IDL: compiled stores unnormalized covariance per trace
    cc_unnormalized = cc_mean_xy .- (mean1 * mean2)
    ac1_unnormalized = ac1_mean_xx .- (mean1^2)
    ac2_unnormalized = ac2_mean_yy .- (mean2^2)
    
    # Compute normalized covariances (CC_ON format) using global means
    # IDL formula (line 311): (total(compiled) / (global_direct_ave * global_delayed_ave)) / n_traces
    # This is: (sum of unnormalized_covariances) / (global_mean1 * global_mean2) / n_traces
    # = (mean of unnormalized_covariances) / (global_mean1 * global_mean2)
    # = (E[xy] - E[x]E[y]) / (E[x]E[y])
    cc_normalized = cc_unnormalized ./ (mean1 * mean2)
    ac1_normalized = ac1_unnormalized ./ (mean1^2)
    ac2_normalized = ac2_unnormalized ./ (mean2^2)
    
    # Debug: Print detailed diagnostic info
    lag0_idx = findfirst(==(0), lags)
    if !isnothing(lag0_idx)
        @info "Cross-correlation diagnostics at lag 0" begin
            println("  mean1 (enhancer ON prob): ", mean1)
            println("  mean2 (gene ON prob): ", mean2)
            println("  E[xy] at lag 0: ", cc_mean_xy[lag0_idx])
            println("  E[x]E[y]: ", mean1 * mean2)
            println("  Unnormalized cov (E[xy] - E[x]E[y]): ", cc_unnormalized[lag0_idx])
            println("  Normalized cov: ", cc_normalized[lag0_idx])
            println("  Number of pairs at lag 0: ", cc_n_pairs[lag0_idx])
            println("  Number of traces: ", length(binary_traces))
        end
        @info "First few lags" begin
            println("  Lags: ", lags[1:min(5, length(lags))])
            println("  CC normalized: ", cc_normalized[1:min(5, length(cc_normalized))])
            println("  CC unnormalized: ", cc_unnormalized[1:min(5, length(cc_unnormalized))])
            println("  E[xy] values: ", cc_mean_xy[1:min(5, length(cc_mean_xy))])
        end
    end
    
    result = (
        cc_normalized=cc_normalized,
        ac1_normalized=ac1_normalized,
        ac2_normalized=ac2_normalized,
        cc_unnormalized=cc_unnormalized,
        ac1_unnormalized=ac1_unnormalized,
        ac2_unnormalized=ac2_unnormalized,
        mean1=mean1,
        mean2=mean2,
        lags=lags
    )
    
    if !bootstrap
        return result
    end
    
    # Bootstrap resampling for confidence intervals
    nthreads = Threads.nthreads()
    cc_norm_bs_threads = [[Float64[] for _ in 1:length(lags)] for _ in 1:nthreads]
    ac1_norm_bs_threads = [[Float64[] for _ in 1:length(lags)] for _ in 1:nthreads]
    ac2_norm_bs_threads = [[Float64[] for _ in 1:length(lags)] for _ in 1:nthreads]
    
    Threads.@threads for b in 1:n_bootstrap
        thread_id = Threads.threadid()
        # Resample traces with replacement (maintains pairing)
        bootstrap_traces = StatsBase.sample(1:n_traces, n_traces, replace=true)
        
        # Compute means for bootstrap sample
        mean1_bs = 0.0
        mean2_bs = 0.0
        n_points_bs = 0
        cc_raw_bs = zeros(length(lags))
        ac1_raw_bs = zeros(length(lags))
        ac2_raw_bs = zeros(length(lags))
        
        for idx in bootstrap_traces
            t = binary_traces[idx]
            cc_raw_bs .+= StatsBase.crosscov(t[:, 1], t[:, 2], lags, demean=false)
            ac1_raw_bs .+= StatsBase.autocov(t[:, 1], lags, demean=false)
            ac2_raw_bs .+= StatsBase.autocov(t[:, 2], lags, demean=false)
            mean1_bs += sum(t[:, 1])
            mean2_bs += sum(t[:, 2])
            n_points_bs += size(t, 1)
        end
        
        mean1_bs /= n_points_bs
        mean2_bs /= n_points_bs
        cc_raw_mean_bs = cc_raw_bs / n_traces
        ac1_raw_mean_bs = ac1_raw_bs / n_traces
        ac2_raw_mean_bs = ac2_raw_bs / n_traces
        
        # Compute normalized covariances for bootstrap sample
        cc_unnorm_bs = cc_raw_mean_bs .- (mean1_bs * mean2_bs)
        ac1_unnorm_bs = ac1_raw_mean_bs .- (mean1_bs^2)
        ac2_unnorm_bs = ac2_raw_mean_bs .- (mean2_bs^2)
        
        cc_norm_bs = cc_unnorm_bs ./ (mean1_bs * mean2_bs)
        ac1_norm_bs = ac1_unnorm_bs ./ (mean1_bs^2)
        ac2_norm_bs = ac2_unnorm_bs ./ (mean2_bs^2)
        
        # Store in thread-local arrays
        for j in 1:length(lags)
            push!(cc_norm_bs_threads[thread_id][j], cc_norm_bs[j])
            push!(ac1_norm_bs_threads[thread_id][j], ac1_norm_bs[j])
            push!(ac2_norm_bs_threads[thread_id][j], ac2_norm_bs[j])
        end
    end
    
    # Combine thread-local arrays
    cc_norm_bs = [Float64[] for _ in 1:length(lags)]
    ac1_norm_bs = [Float64[] for _ in 1:length(lags)]
    ac2_norm_bs = [Float64[] for _ in 1:length(lags)]
    for t in 1:nthreads
        for j in 1:length(lags)
            append!(cc_norm_bs[j], cc_norm_bs_threads[t][j])
            append!(ac1_norm_bs[j], ac1_norm_bs_threads[t][j])
            append!(ac2_norm_bs[j], ac2_norm_bs_threads[t][j])
        end
    end
    
    # Compute percentiles
    cc_norm_lower = [quantile(cc_norm_bs[j], 0.025) for j in 1:length(lags)]
    cc_norm_median = [median(cc_norm_bs[j]) for j in 1:length(lags)]
    cc_norm_upper = [quantile(cc_norm_bs[j], 0.975) for j in 1:length(lags)]
    cc_norm_std = [std(cc_norm_bs[j]) for j in 1:length(lags)]
    
    ac1_norm_lower = [quantile(ac1_norm_bs[j], 0.025) for j in 1:length(lags)]
    ac1_norm_median = [median(ac1_norm_bs[j]) for j in 1:length(lags)]
    ac1_norm_upper = [quantile(ac1_norm_bs[j], 0.975) for j in 1:length(lags)]
    ac1_norm_std = [std(ac1_norm_bs[j]) for j in 1:length(lags)]
    
    ac2_norm_lower = [quantile(ac2_norm_bs[j], 0.025) for j in 1:length(lags)]
    ac2_norm_median = [median(ac2_norm_bs[j]) for j in 1:length(lags)]
    ac2_norm_upper = [quantile(ac2_norm_bs[j], 0.975) for j in 1:length(lags)]
    ac2_norm_std = [std(ac2_norm_bs[j]) for j in 1:length(lags)]
    
    return (
        cc_normalized=cc_normalized,
        ac1_normalized=ac1_normalized,
        ac2_normalized=ac2_normalized,
        cc_unnormalized=cc_unnormalized,
        ac1_unnormalized=ac1_unnormalized,
        ac2_unnormalized=ac2_unnormalized,
        cc_normalized_lower=cc_norm_lower,
        cc_normalized_median=cc_norm_median,
        cc_normalized_upper=cc_norm_upper,
        cc_normalized_std=cc_norm_std,
        ac1_normalized_lower=ac1_norm_lower,
        ac1_normalized_median=ac1_norm_median,
        ac1_normalized_upper=ac1_norm_upper,
        ac1_normalized_std=ac1_norm_std,
        ac2_normalized_lower=ac2_norm_lower,
        ac2_normalized_median=ac2_norm_median,
        ac2_normalized_upper=ac2_norm_upper,
        ac2_normalized_std=ac2_norm_std,
        mean1=mean1,
        mean2=mean2,
        lags=lags
    )
end

"""
    score_model_on_state(unit1_filename, unit2_filename, r, transitions, G, R, S, insertstep, coupling, interval, probfn, lags; trace_id_col=nothing, n_bootstrap=10000)

Score model predictions against data for ON state cross-covariances and autocovariances.

This function walks through the complete stack:
1. Loads data from two CSV files (unit1 and unit2) using `load_predicted_traces_csv`
2. Computes ON state correlations from data using `compute_binary_correlations`
3. Gets theoretical predictions from `covariance_functions`
4. Compares them and computes scoring metrics

# Arguments
- `unit1_filename::String`: Path to CSV file for unit 1 (enhancer)
- `unit2_filename::String`: Path to CSV file for unit 2 (gene)
- `r`: Rate parameters
- `transitions`: Transition definitions
- `G`: Gene definitions
- `R`: Rate definitions
- `S`: State definitions
- `insertstep`: Insert step definitions
- `coupling`: Coupling parameters
- `interval`: Time interval
- `probfn`: Probability function
- `lags::Vector{Int}`: Time lags for correlation calculation
- `trace_id_col::Union{String, Nothing}=nothing`: Column name for trace IDs (if traces are grouped)
- `n_bootstrap::Int=10000`: Number of bootstrap replicates for confidence intervals

# Returns
Named tuple with:
- `cc_theory`: Theoretical normalized ON state cross-covariance
- `cc_data`: Empirical normalized ON state cross-covariance
- `cc_data_lower`, `cc_data_upper`: Bootstrap confidence intervals for data
- `ac1_theory`, `ac2_theory`: Theoretical normalized ON state autocovariances
- `ac1_data`, `ac2_data`: Empirical normalized ON state autocovariances
- `ac1_data_lower`, `ac1_data_upper`: Bootstrap CIs for AC1
- `ac2_data_lower`, `ac2_data_upper`: Bootstrap CIs for AC2
- `cc_l2_norm`, `cc_linf_norm`: L² and L∞ norms for cross-covariance
- `ac1_l2_norm`, `ac1_linf_norm`: L² and L∞ norms for AC1
- `ac2_l2_norm`, `ac2_linf_norm`: L² and L∞ norms for AC2
- `cc_z_scores`: Z-scores for cross-covariance (theory - data) / data_std
- `cc_coverage`: Boolean array indicating if theory is within bootstrap CI
- `cc_coverage_fraction`: Fraction of lags where theory is within CI
- `mean1_data`, `mean2_data`: Empirical ON state means (fraction ON)
- `mean1_theory`, `mean2_theory`: Theoretical ON state means
- `lags`: Time lags

# Example
```julia
results = score_model_on_state(
    "predictedtraces_trace-HBEC-nstate-h_enhancer_MYC_3301_1.csv",
    "predictedtraces_trace-HBEC-nstate-h_gene_MYC_3301_1.csv",
    r, transitions, G, R, S, insertstep, coupling, interval, probfn, lags
)
```
"""
function score_model_on_state(unit1_filename, unit2_filename, r, transitions, G, R, S, insertstep, coupling, interval, probfn, lags; trace_id_col=nothing, n_bootstrap=10000)
    # Step 1: Load data from CSV files
    data_traces = load_predicted_traces_csv(unit1_filename, unit2_filename; trace_id_col=trace_id_col)
    
    # Step 2: Compute ON state correlations from data (with bootstrap)
    data_result = compute_binary_correlations(data_traces, lags; bootstrap=true, n_bootstrap=n_bootstrap)
    
    # Step 3: Get theoretical predictions
    # covariance_functions expects positive lags only, but returns full symmetric array
    positive_lags = unique([abs.(lags); lags[lags .>= 0]])
    sort!(positive_lags)
    
    ac1_theory, ac2_theory, cc_intensity_theory, ccON_theory, full_lags, m1, m2, v1, v2, mON1, mON2, ac1ON_theory, ac2ON_theory = 
        StochasticGene.covariance_functions(r, transitions, G, R, S, insertstep, interval, probfn, coupling, positive_lags)
    
    # Match lags: theory returns full symmetric array, we need to extract matching lags
    # Theory returns: [negative_lags..., 0, positive_lags...]
    # We need to find indices for our requested lags
    ccON_theory_matched = zeros(length(lags))
    ac1ON_theory_matched = zeros(length(lags))
    ac2ON_theory_matched = zeros(length(lags))
    
    for (i, lag) in enumerate(lags)
        idx = findfirst(==(lag), full_lags)
        if !isnothing(idx)
            ccON_theory_matched[i] = ccON_theory[idx]
            ac1ON_theory_matched[i] = ac1ON_theory[idx]
            ac2ON_theory_matched[i] = ac2ON_theory[idx]
        else
            # Linear interpolation if lag not found
            idx1 = findlast(<(lag), full_lags)
            idx2 = findfirst(>(lag), full_lags)
            if !isnothing(idx1) && !isnothing(idx2)
                w1 = (full_lags[idx2] - lag) / (full_lags[idx2] - full_lags[idx1])
                w2 = (lag - full_lags[idx1]) / (full_lags[idx2] - full_lags[idx1])
                ccON_theory_matched[i] = w1 * ccON_theory[idx1] + w2 * ccON_theory[idx2]
                ac1ON_theory_matched[i] = w1 * ac1ON_theory[idx1] + w2 * ac1ON_theory[idx2]
                ac2ON_theory_matched[i] = w1 * ac2ON_theory[idx1] + w2 * ac2ON_theory[idx2]
            end
        end
    end
    
    # Step 4: Compute scoring metrics
    # Cross-covariance metrics
    cc_l2_norm = sqrt(sum((ccON_theory_matched .- data_result.cc_normalized).^2))
    cc_linf_norm = maximum(abs.(ccON_theory_matched .- data_result.cc_normalized))
    
    # Z-scores: (theory - data_mean) / data_std
    cc_z_scores = (ccON_theory_matched .- data_result.cc_normalized) ./ data_result.cc_normalized_se
    cc_z_scores[isnan.(cc_z_scores)] .= 0.0  # Handle division by zero
    
    # Coverage: whether theory is within bootstrap CI
    cc_coverage = (ccON_theory_matched .>= data_result.cc_normalized_lower) .& 
                  (ccON_theory_matched .<= data_result.cc_normalized_upper)
    
    # Autocovariance metrics
    ac1_l2_norm = sqrt(sum((ac1ON_theory_matched .- data_result.ac1_normalized).^2))
    ac1_linf_norm = maximum(abs.(ac1ON_theory_matched .- data_result.ac1_normalized))
    
    ac2_l2_norm = sqrt(sum((ac2ON_theory_matched .- data_result.ac2_normalized).^2))
    ac2_linf_norm = maximum(abs.(ac2ON_theory_matched .- data_result.ac2_normalized))
    
    return (
        # Cross-covariance
        cc_theory=ccON_theory_matched,
        cc_data=data_result.cc_normalized,
        cc_data_lower=data_result.cc_normalized_lower,
        cc_data_median=data_result.cc_normalized_median,
        cc_data_upper=data_result.cc_normalized_upper,
        cc_data_se=data_result.cc_normalized_se,
        cc_l2_norm=cc_l2_norm,
        cc_linf_norm=cc_linf_norm,
        cc_z_scores=cc_z_scores,
        cc_coverage=cc_coverage,
        cc_coverage_fraction=mean(cc_coverage),
        # Autocovariance 1
        ac1_theory=ac1ON_theory_matched,
        ac1_data=data_result.ac1_normalized,
        ac1_data_lower=data_result.ac1_normalized_lower,
        ac1_data_median=data_result.ac1_normalized_median,
        ac1_data_upper=data_result.ac1_normalized_upper,
        ac1_data_se=data_result.ac1_normalized_se,
        ac1_l2_norm=ac1_l2_norm,
        ac1_linf_norm=ac1_linf_norm,
        # Autocovariance 2
        ac2_theory=ac2ON_theory_matched,
        ac2_data=data_result.ac2_normalized,
        ac2_data_lower=data_result.ac2_normalized_lower,
        ac2_data_median=data_result.ac2_normalized_median,
        ac2_data_upper=data_result.ac2_normalized_upper,
        ac2_data_se=data_result.ac2_normalized_se,
        ac2_l2_norm=ac2_l2_norm,
        ac2_linf_norm=ac2_linf_norm,
        # Means
        mean1_data=data_result.mean1,
        mean2_data=data_result.mean2,
        mean1_theory=mON1,
        mean2_theory=mON2,
        lags=lags
    )
end

"""
    score_models_from_crosscovariance_files(folder_path, r, transitions, G, R, S, insertstep, interval, probfn, lags; pattern="crosscovariance_tracejoint-HBEC-nstate_enhancer-gene")

Score model predictions against data from precomputed crosscovariance CSV files.

This function:
1. Reads all crosscovariance CSV files from a folder matching the pattern
2. Extracts coupling model identifier (e.g., "31", "R2") from filenames
3. Reads `tau` and `cc_ON` columns from each file
4. Computes theoretical predictions for each coupling model
5. Scores each model against its corresponding data

# Arguments
- `folder_path::String`: Path to folder containing crosscovariance CSV files
- `r`: Rate parameters (can be a function of coupling model or a single value)
- `transitions`: Transition definitions
- `G`: Gene definitions
- `R`: Rate definitions
- `S`: State definitions
- `insertstep`: Insert step definitions
- `interval`: Time interval
- `probfn`: Probability function
- `coupling`: Coupling parameter (can be a single value, a Dict mapping coupling_model_str to coupling, or a Function)
- `lags::Vector{Int}`: Time lags for correlation calculation (should match tau in files)
- `pattern::String="crosscovariance_tracejoint-HBEC-nstate_enhancer-gene"`: Filename pattern to match

# Returns
Dictionary mapping coupling model identifiers to named tuples with:
- `coupling_model::String`: Coupling model identifier
- `cc_theory`: Theoretical normalized ON state cross-covariance
- `cc_data`: Empirical normalized ON state cross-covariance from file
- `tau`: Time lags from file
- `cc_l2_norm`: L² norm (goodness of fit)
- `cc_linf_norm`: L∞ norm (worst-case error)
- `cc_z_scores`: Z-scores (if theory std available)
- `mean1_theory`, `mean2_theory`: Theoretical ON state means

# Example
```julia
# If coupling is the same for all models:
results = score_models_from_crosscovariance_files(
    "/path/to/folder",
    r, transitions, G, R, S, insertstep, interval, probfn, coupling, lags
)

# If coupling varies by model (as Dict):
coupling_dict = Dict("31" => coupling31, "R2" => couplingR2)
results = score_models_from_crosscovariance_files(
    "/path/to/folder",
    r, transitions, G, R, S, insertstep, interval, probfn, coupling_dict, lags
)

# Access results for coupling model "31":
results["31"].cc_l2_norm
```
"""
function score_models_from_crosscovariance_files(folder_path, r, transitions, G, R, S, insertstep, interval, probfn, coupling, lags; pattern="crosscovariance_tracejoint-HBEC-nstate_enhancer-gene")
    # Find all matching files
    files = filter(f -> startswith(f, pattern) && endswith(f, ".csv"), readdir(folder_path))
    
    if isempty(files)
        error("No files matching pattern '$pattern*.csv' found in folder '$folder_path'")
    end
    
    results = Dict{String, NamedTuple}()
    
    for file in files
        # Extract coupling model from filename
        # Pattern: crosscovariance_tracejoint-HBEC-nstate_enhancer-gene{COUPLING}_MYC_...
        # Find position after "gene"
        gene_pos = findfirst("gene", file)
        if isnothing(gene_pos)
            @warn "Could not find 'gene' in filename '$file', skipping"
            continue
        end
        
        # Extract two characters after "gene"
        coupling_start = gene_pos.stop + 1
        if coupling_start + 1 > length(file)
            @warn "Filename '$file' too short to extract coupling model, skipping"
            continue
        end
        
        coupling_model_str = file[coupling_start:coupling_start+1]
        
        # Read CSV file (CSV and DataFrames should be available in the module)
        filepath = joinpath(folder_path, file)
        df = CSV.read(filepath, DataFrame)
        
        # Check required columns
        if !("tau" in names(df)) || !("cc_ON" in names(df))
            @warn "File '$file' missing required columns 'tau' or 'cc_ON', skipping"
            continue
        end
        
        # Extract data
        tau_data = df.tau
        cc_ON_data = df.cc_ON
        
        # Match lags: we need to find which lags from our requested lags match tau_data
        # Interpolate if needed
        cc_ON_matched = zeros(length(lags))
        for (i, lag) in enumerate(lags)
            idx = findfirst(==(lag), tau_data)
            if !isnothing(idx)
                cc_ON_matched[i] = cc_ON_data[idx]
            else
                # Linear interpolation if lag not found
                idx1 = findlast(<(lag), tau_data)
                idx2 = findfirst(>(lag), tau_data)
                if !isnothing(idx1) && !isnothing(idx2)
                    w1 = (tau_data[idx2] - lag) / (tau_data[idx2] - tau_data[idx1])
                    w2 = (lag - tau_data[idx1]) / (tau_data[idx2] - tau_data[idx1])
                    cc_ON_matched[i] = w1 * cc_ON_data[idx1] + w2 * cc_ON_data[idx2]
                else
                    @warn "Lag $lag not found in file '$file' and cannot interpolate, using 0"
                end
            end
        end
        
        # Get theoretical predictions
        # Handle coupling parameter: can be single value, Dict, or Function
        coupling_actual = if isa(coupling, Dict)
            get(coupling, coupling_model_str, error("Coupling model '$coupling_model_str' not found in coupling dictionary"))
        elseif isa(coupling, Function)
            coupling(coupling_model_str)
        else
            coupling  # Single value, use for all models
        end
        
        # Handle r parameter: can be single value or Function
        r_actual = isa(r, Function) ? r(coupling_model_str) : r
        
        positive_lags = unique([abs.(lags); lags[lags .>= 0]])
        sort!(positive_lags)
        
        ac1_theory, ac2_theory, cc_intensity_theory, ccON_theory, full_lags, m1, m2, v1, v2, mON1, mON2, ac1ON_theory, ac2ON_theory = 
            StochasticGene.covariance_functions(r_actual, transitions, G, R, S, insertstep, interval, probfn, coupling_actual, positive_lags)
        
        # Match theoretical lags to requested lags
        ccON_theory_matched = zeros(length(lags))
        for (i, lag) in enumerate(lags)
            idx = findfirst(==(lag), full_lags)
            if !isnothing(idx)
                ccON_theory_matched[i] = ccON_theory[idx]
            else
                # Linear interpolation
                idx1 = findlast(<(lag), full_lags)
                idx2 = findfirst(>(lag), full_lags)
                if !isnothing(idx1) && !isnothing(idx2)
                    w1 = (full_lags[idx2] - lag) / (full_lags[idx2] - full_lags[idx1])
                    w2 = (lag - full_lags[idx1]) / (full_lags[idx2] - full_lags[idx1])
                    ccON_theory_matched[i] = w1 * ccON_theory[idx1] + w2 * ccON_theory[idx2]
                end
            end
        end
        
        # Compute scoring metrics
        cc_l2_norm = sqrt(sum((ccON_theory_matched .- cc_ON_matched).^2))
        cc_linf_norm = maximum(abs.(ccON_theory_matched .- cc_ON_matched))
        
        # Store results
        results[coupling_model_str] = (
            coupling_model=coupling_model_str,
            cc_theory=ccON_theory_matched,
            cc_data=cc_ON_matched,
            tau=tau_data,
            cc_l2_norm=cc_l2_norm,
            cc_linf_norm=cc_linf_norm,
            mean1_theory=mON1,
            mean2_theory=mON2,
            lags=lags
        )
    end
    
    return results
end

"""
    score_models_from_traces(enhancer_file, gene_file, crosscov_folder; 
        crosscov_pattern="crosscovariance_tracejoint-HBEC-nstate_enhancer-gene")

Score model predictions against empirical cross-covariance computed from trace predictions.

The crosscovariance files contain the theoretical predictions (precomputed).
Workflow:
1. Load enhancer and gene trace prediction files (from separate fits)
2. Read theoretical cross-covariance from CSV files in crosscov_folder (includes tau/lags)
3. Compute empirical cross-covariance (Xcor) from traces using lags from the file
4. Score theoretical predictions (from files) against empirical data (from traces)

The coupling model is extracted from crosscovariance filenames (two characters after "gene").
Lags are read from the `tau` column in the crosscovariance files.

# Arguments
- `enhancer_file::String`: Path to enhancer prediction CSV file
- `gene_file::String`: Path to gene prediction CSV file
- `crosscov_folder::String`: Folder containing crosscovariance CSV files (with theoretical predictions)
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
function score_models_from_traces(enhancer_file, gene_file, crosscov_folder;
    crosscov_pattern="crosscovariance_tracejoint-HBEC-nstate_enhancer-gene")
    
    # Find all crosscovariance files (these contain the coupling model identifier)
    crosscov_files = filter(f -> startswith(f, crosscov_pattern) && endswith(f, ".csv"), readdir(crosscov_folder))
    
    if isempty(crosscov_files)
        error("No crosscovariance files matching pattern '$crosscov_pattern*.csv' found in folder '$crosscov_folder'")
    end
    
    results = Dict{String, NamedTuple}()
    
    # 1. Load traces first to determine valid lag range
    data_traces = try
        load_predicted_traces_csv(enhancer_file, gene_file)
    catch e
        error("Error loading traces from '$enhancer_file' and '$gene_file': $e")
    end
    
    # Determine minimum trace length
    min_trace_length = minimum([size(t, 1) for t in data_traces])
    
    # Get lags from the first crosscovariance file (assume all files use same lags)
    first_file = joinpath(crosscov_folder, crosscov_files[1])
    df_first = CSV.read(first_file, DataFrame)
    if !("tau" in names(df_first))
        error("First crosscovariance file missing 'tau' column")
    end
    lags_all = Vector{Int}(df_first.tau)
    
    # Filter lags to only include those valid for trace length
    # Use maximum possible lags: abs(lag) < min_trace_length (strictly less than)
    # This maximizes the lag range in both directions
    max_valid_lag = min_trace_length - 1
    lags_all_filtered = filter(lag -> abs(lag) < max_valid_lag, lags_all)
    
    # Optionally use higher resolution: if lags are coarsely spaced, interpolate
    # Check if we should use higher resolution (lag step > 1)
    if length(lags_all_filtered) > 1
        lag_step = abs(lags_all_filtered[2] - lags_all_filtered[1])
        if lag_step > 1
            # Create higher resolution lags
            min_lag = minimum(lags_all_filtered)
            max_lag = maximum(lags_all_filtered)
            lags = collect(min_lag:max_lag)  # Step size of 1
            @info "Using higher resolution: original step=$lag_step, new step=1, lags from $min_lag to $max_lag"
        else
            lags = lags_all_filtered
        end
    else
        lags = lags_all_filtered
    end
    
    if isempty(lags)
        error("No valid lags found. Trace length is $min_trace_length, but lags range from $(minimum(lags_all)) to $(maximum(lags_all))")
    end
    
    # Compute empirical correlations using filtered lags
    data_result = try
        compute_binary_correlations(data_traces, lags; bootstrap=false)
    catch e
        error("Error computing correlations from traces: $e")
    end
    
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
        crosscov_filepath = joinpath(crosscov_folder, crosscov_file)
        
        # 2. Read crosscovariance from file
        try
            df = CSV.read(crosscov_filepath, DataFrame)
            if !("tau" in names(df)) || !("cc_ON" in names(df))
                @warn "File '$crosscov_file' missing required columns 'tau' or 'cc_ON', skipping"
                continue
            end
            
            tau_data = Vector{Int}(df.tau)
            cc_ON_theory_all = df.cc_ON  # Theoretical predictions from file
            
            # Filter theory to only valid lags (matching what we used for empirical)
            cc_ON_theory = Float64[]
            for lag in lags
                idx = findfirst(==(lag), tau_data)
                if !isnothing(idx)
                    push!(cc_ON_theory, cc_ON_theory_all[idx])
                else
                    # Linear interpolation if lag not found
                    idx1 = findlast(<(lag), tau_data)
                    idx2 = findfirst(>(lag), tau_data)
                    if !isnothing(idx1) && !isnothing(idx2)
                        w1 = (tau_data[idx2] - lag) / (tau_data[idx2] - tau_data[idx1])
                        w2 = (lag - tau_data[idx1]) / (tau_data[idx2] - tau_data[idx1])
                        val = w1 * cc_ON_theory_all[idx1] + w2 * cc_ON_theory_all[idx2]
                        push!(cc_ON_theory, val)
                    else
                        @warn "Could not find or interpolate lag $lag in file '$crosscov_file', skipping"
                    end
                end
            end
            
            if length(cc_ON_theory) != length(lags)
                @warn "Mismatch in filtered theory length for '$crosscov_file', skipping"
                continue
            end
            
            cc_ON_theory = Vector{Float64}(cc_ON_theory)
            
            # 3. Compute scoring metrics (theory from file vs empirical from traces)
            cc_l2_norm = sqrt(sum((cc_ON_theory .- data_result.cc_normalized).^2))
            cc_linf_norm = maximum(abs.(cc_ON_theory .- data_result.cc_normalized))
            
            # Store results
            results[coupling_model_str] = (
                coupling_model=coupling_model_str,
                cc_empirical=data_result.cc_normalized,
                cc_theory=cc_ON_theory,
                cc_l2_norm=cc_l2_norm,
                cc_linf_norm=cc_linf_norm,
                mean1_empirical=data_result.mean1,
                mean2_empirical=data_result.mean2,
                lags=lags
            )
            
        catch e
            @warn "Error reading crosscovariance file '$crosscov_file': $e, skipping"
            continue
        end
    end
    
    return results
end

"""
    summarize_model_scores(results::Dict{String, NamedTuple})

Print a summary table of model scores sorted by L² norm (best fit first).
"""
function summarize_model_scores(results::Dict{String, NamedTuple})
    println("\nModel Scoring Summary (sorted by L² norm, best fit first):")
    println("=" ^ 80)
    println(@sprintf("%-10s %15s %15s %15s", "Model", "L² Norm", "L∞ Norm", "Mean1"))
    println("-" ^ 80)
    
    # Sort by L² norm
    sorted_models = sort(collect(keys(results)), by=k -> results[k].cc_l2_norm)
    
    for model in sorted_models
        r = results[model]
        println(@sprintf("%-10s %15.6f %15.6f %15.6f", 
            model, r.cc_l2_norm, r.cc_linf_norm, r.mean1_empirical))
    end
    println("=" ^ 80)
    println("\nNote: cc_empirical is the same for all models (computed from data).")
    println("      Differences are in cc_theory, which determines the norms.\n")
end

