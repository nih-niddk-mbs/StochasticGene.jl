# load_data Function

Load experimental data from files into StochasticGene.jl data structures.

## Syntax

```julia
load_data(datatype, dttype, datapath, label, gene, datacond, traceinfo, temprna, datacol=3, zeromedian=false)
```

## Arguments

### Required Arguments

- `datatype::String`: Type of data to load
- `dttype::Vector`: Dwell time types (for dwell time data)
- `datapath::String`: Path to data file(s)
- `label::String`: Label for the dataset
- `gene::String`: Gene name
- `datacond::String`: Experimental condition
- `traceinfo::Tuple`: Trace parameters (for trace data)
- `temprna::Int`: Temperature scaling factor for RNA data

### Optional Arguments

- `datacol::Int = 3`: Column index for trace data
- `zeromedian::Bool = false`: Whether to zero-center traces

## Supported Data Types

### RNA Data Types

#### "rna" - RNA Count Histograms
Load steady-state RNA count distributions from smFISH or scRNA-seq data.

**File Format**: Text files with RNA count histograms
**Returns**: `RNAData` structure

```julia
# Load RNA histogram data
data = load_data(
    "rna",                    # datatype
    String[],                 # dttype (unused for RNA)
    "data/smFISH/",          # datapath
    "control",               # label
    "MYC",                   # gene
    "MOCK",                  # datacond
    (),                      # traceinfo (unused for RNA)
    1                        # temprna
)
```

#### "rnacount" - Individual RNA Counts
Load individual RNA count measurements with yield factors.

**File Format**: Text files with individual counts and yield factors
**Returns**: `RNACountData` structure

```julia
# Load individual RNA counts
data = load_data(
    "rnacount",
    String[],
    "data/single_cell/",
    "experiment1",
    "ACTB",
    "treated",
    (),
    1
)
```

#### "rnaonoff" - RNA + ON/OFF Times
Load combined RNA counts and ON/OFF state duration data.

**File Format**: Two files - RNA histogram and ON/OFF times
**Returns**: `RNAOnOffData` structure

```julia
# Load RNA + ON/OFF data
data = load_data(
    "rnaonoff",
    String[],
    ["data/rna_hist.txt", "data/onoff_times.txt"],  # Multiple files
    "live_cell",
    "SOX2",
    "control",
    (),
    2  # Divide histogram by 2
)
```

#### "rnadwelltime" - RNA + Dwell Times
Load RNA counts with dwell time distributions.

**File Format**: Multiple files - RNA histogram and dwell time data
**Returns**: `RNADwellTimeData` structure

```julia
# Load RNA + dwell time data
data = load_data(
    "rnadwelltime",
    ["ON", "OFF"],           # Dwell time types
    ["data/rna.txt", "data/dwellON.txt", "data/dwellOFF.txt"],
    "combined",
    "GENE1",
    "condition1",
    (),
    1
)
```

### Trace Data Types

#### "trace" - Intensity Traces
Load fluorescence intensity time series data.

**File Format**: Text files with intensity traces
**Returns**: `TraceData` structure

```julia
# Load intensity traces
data = load_data(
    "trace",
    String[],
    "data/traces/",
    "live_imaging",
    "MYC",
    "control",
    (1.0, 1.0, -1, 1.0),    # traceinfo: (interval, start, end, weight)
    1,
    3,                       # datacol: use column 3
    false                    # zeromedian: don't zero-center
)
```

#### "tracerna" - Traces + RNA Histogram
Load both intensity traces and RNA histogram data.

**File Format**: Two directories - trace files and RNA histogram
**Returns**: `TraceRNAData` structure

```julia
# Load traces + RNA data
data = load_data(
    "tracerna",
    String[],
    ["data/traces/", "data/rna/"],  # Two data paths
    "combined_exp",
    "GENE2",
    "treated",
    (0.5, 0.0, 100.0, 0.8),        # 0.5 min intervals, 0-100 min, 80% weight
    1,
    3,
    true                            # Zero-center traces
)
```

#### "tracejoint" - Simultaneous Traces
Load traces recorded simultaneously (e.g., different fluorescent channels).

**File Format**: Multi-column trace files
**Returns**: `TraceData` structure

```julia
# Load joint traces
data = load_data(
    "tracejoint",
    String[],
    "data/joint_traces/",
    "dual_channel",
    "GENE_TF",
    "stimulus",
    (2.0, 0.0, -1, [0.9, 0.85]),   # Different weights per channel
    1,
    [3, 4],                        # Use columns 3 and 4
    [true, false]                  # Zero-center first channel only
)
```

#### "tracegrid" - Grid-Based Traces
Load traces for grid-based parameter exploration.

**File Format**: Structured trace files
**Returns**: `TraceData` structure

```julia
# Load grid traces
data = load_data(
    "tracegrid",
    String[],
    "data/grid_traces/",
    "parameter_grid",
    "TARGET",
    "sweep",
    (1.0, 0.0, 50.0, 1.0),
    1
)
```

### Dwell Time Data

#### "dwelltime" - Dwell Time Distributions Only
Load dwell time distributions without RNA data.

**File Format**: Text files with dwell time histograms
**Returns**: `DwellTimeData` structure

```julia
# Load dwell time data
data = load_data(
    "dwelltime",
    ["ON", "OFF", "BURST"],
    ["data/dwell_ON.txt", "data/dwell_OFF.txt", "data/dwell_BURST.txt"],
    "dwell_analysis",
    "GENE3",
    "control",
    (),
    1
)
```

## Trace Parameters (traceinfo)

For trace data types, the `traceinfo` tuple contains:

1. **interval** (`Float64`): Time interval between frames (minutes)
2. **start** (`Float64`): Start time or frame number
3. **end** (`Float64`): End time or frame number (-1 for last frame)
4. **weight** (`Float64` or `Vector{Float64}`): Fraction of active traces or per-channel weights
5. **background** (`Vector` or `Vector{Vector}`, optional): Background values per trace

```julia
# Basic trace info
traceinfo = (1.0, 1.0, -1, 1.0)  # 1 min intervals, frame 1 to end, 100% weight

# Multi-channel with backgrounds
traceinfo = (
    0.5,                    # 30 second intervals
    0.0,                    # Start at time 0
    100.0,                  # End at 100 minutes
    [0.8, 0.9],            # 80% weight channel 1, 90% channel 2
    [[10.0, 15.0], [5.0, 8.0]]  # Background per channel per trace
)
```

## File Formats

### RNA Histogram Files
```
# RNA count histogram
# Count  Frequency
0       150
1       200
2       180
3       120
4       80
...
```

### Trace Files
```
# Intensity trace data
# Time   Background   Intensity   State
0.0     10.0         25.3        1
1.0     10.2         28.1        1
2.0     9.8          24.5        1
3.0     10.1         45.2        2
...
```

### Dwell Time Files
```
# Dwell time histogram
# Time_bin   Count
0.5         25
1.0         40
1.5         35
2.0         30
...
```

## Examples

### Complete RNA Analysis Workflow

```julia
using StochasticGene

# Load RNA data
rna_data = load_data(
    "rna",
    String[],
    "data/HCT116_smFISH/",
    "HCT116_control",
    "MYC",
    "MOCK",
    (),
    1
)

# Display data properties
println("Gene: ", rna_data.gene)
println("RNA range: ", length(rna_data.histRNA))
println("Total counts: ", sum(rna_data.histRNA))
```

### Multi-Modal Data Loading

```julia
# Load combined RNA and trace data
combined_data = load_data(
    "tracerna",
    String[],
    ["data/traces/MYC_traces/", "data/rna/MYC_histogram.txt"],
    "multimodal",
    "MYC",
    "control",
    (1.0, 0.0, 120.0, 0.85),
    1,
    3,
    true
)

# Access both data types
traces = combined_data.trace
rna_hist = combined_data.histRNA
```

### Hierarchical Data Loading

```julia
# Load multiple conditions for hierarchical analysis
conditions = ["control", "treated1", "treated2"]
all_data = []

for cond in conditions
    data = load_data(
        "trace",
        String[],
        "data/traces/",
        "hierarchical_$cond",
        "SOX2",
        cond,
        (0.5, 0.0, 60.0, 0.9),
        1,
        3,
        true
    )
    push!(all_data, data)
end
```

### Custom Data Processing

```julia
# Load with custom processing
function custom_load_processing(datapath, gene, condition)
    # Custom preprocessing logic
    raw_data = load_data(
        "trace",
        String[],
        datapath,
        "custom_$condition",
        gene,
        condition,
        (2.0, 5.0, 100.0, 0.95),
        1,
        3,
        true
    )
    
    # Additional processing
    processed_traces = []
    for trace in raw_data.trace[1]
        # Apply custom filtering, normalization, etc.
        filtered_trace = smooth_trace(trace)
        push!(processed_traces, filtered_trace)
    end
    
    return raw_data, processed_traces
end
```

## Error Handling

The function includes comprehensive error checking:

```julia
# Invalid data type
try
    data = load_data("invalid_type", ...)
catch ArgumentError as e
    println("Error: ", e.msg)
end

# Missing files
try
    data = load_data("rna", String[], "nonexistent/path/", ...)
catch SystemError as e
    println("File error: ", e.msg)
end

# Inconsistent trace info
try
    data = load_data("trace", String[], path, label, gene, cond, 
                    (1.0, 100.0, 50.0, 1.0), ...)  # end < start
catch ArgumentError as e
    println("Invalid trace parameters: ", e.msg)
end
```

## Performance Considerations

1. **Large Files**: Use appropriate data types for memory efficiency
2. **Multiple Files**: Consider loading in batches for very large datasets
3. **Trace Data**: Zero-centering can reduce memory usage for fitting
4. **File Formats**: CSV files are slower than tab-delimited files

## Data Validation

The function validates:
- File existence and readability
- Data format consistency
- Parameter compatibility
- Trace parameter validity
- Column specifications for trace data

## See Also

- [`RNAData`](@ref): RNA data structure
- [`TraceData`](@ref): Trace data structure
- [`fit`](@ref): Model fitting function
- [`simulator`](@ref): Data simulation function
- [`write_dataframes`](@ref): Write data to files