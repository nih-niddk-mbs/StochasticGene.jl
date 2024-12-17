# This file is part of StochasticGene.jl   

###  data_types.jl

# Data types
"""
    AbstractExperimentalData

Abstract type for experimental data
"""
abstract type AbstractExperimentalData end
"""
    AbstractSampleData

Abstract type for data in the form of samples
"""
abstract type AbstractSampleData <: AbstractExperimentalData end
"""
    AbstractHistogramData

Abstract type for data in the form of a histogram (probability distribution)
"""
abstract type AbstractHistogramData <: AbstractExperimentalData end
"""
    AbstractTraceData
    
Abstract type for intensity time series data
"""
abstract type AbstractTraceData <: AbstractExperimentalData end

"""
  AbstractRNAData{hType}

Abstract type for steady state RNA histogram data
"""
abstract type AbstractRNAData{hType} <: AbstractHistogramData end


"""
    AbstractTraceHistogramData
    
Abstract type for intensity time series data with RNA histogram data
"""
abstract type AbstractTraceHistogramData <: AbstractExperimentalData end

# Data structures 
#
# Do not use underscore "_" in label

"""
    RNAData{nType,hType}

Structure for storing RNA histogram data.

# Fields
- `label::String`: Label for the data set.
- `gene::String`: Gene name (case sensitive).
- `nRNA::nType`: Length of the histogram.
- `histRNA::hType`: RNA histograms.
"""
struct RNAData{nType,hType} <: AbstractRNAData{hType}
    label::String
    gene::String
    nRNA::nType
    histRNA::hType
end

struct DwellTimeData <: AbstractHistogramData
    label::String
    gene::String
    bins::Vector
    DwellTimes::Vector
    DTtypes::Vector
end





"""
    RNAOnOffData

Structure for storing RNA ON/OFF time data.

# Fields
- `label::String`: Label for the data set.
- `gene::String`: Gene name (case sensitive).
- `nRNA::Int`: Length of the histogram.
- `histRNA::Vector`: RNA histograms.
- `bins::Vector`: Number of live cell recording time bins.
- `ON::Vector`: ON time probability density.
- `OFF::Vector`: OFF time probability density.
"""
struct RNAOnOffData <: AbstractHistogramData
    label::String
    gene::String
    nRNA::Int
    histRNA::Vector
    bins::Vector
    ON::Vector
    OFF::Vector
end
"""
    RNADwellTimeData

Structure for storing RNA dwell time data.

# Fields
- `label::String`: Label for the data set.
- `gene::String`: Gene name (case sensitive).
- `nRNA::Int`: Length of the histogram.
- `histRNA::Array`: RNA histograms.
- `bins::Vector{Vector}`: Number of live cell recording time bins.
- `DwellTimes::Vector{Vector}`: Dwell times.
- `DTtypes::Vector`: Types of dwell times.
"""
struct RNADwellTimeData <: AbstractHistogramData
    label::String
    gene::String
    nRNA::Int
    histRNA::Array
    bins::Vector{Vector}
    DwellTimes::Vector{Vector}
    DTtypes::Vector
end
"""
    TraceData{labelType,geneType,traceType}

Structure for storing trace data.

# Fields
- `label::labelType`: Label for the data set.
- `gene::geneType`: Gene name (case sensitive).
- `interval::Float64`: Time interval between trace points.
- `trace::traceType`: Trace data.
"""
struct TraceData{labelType,geneType,traceType} <: AbstractTraceData
    label::labelType
    gene::geneType
    interval::Float64
    trace::traceType
end
"""
    TraceRNAData{traceType,hType}

Structure for storing trace RNA histogram data.

# Fields
- `label::String`: Label for the data set.
"""
struct TraceRNAData{traceType,hType} <: AbstractTraceHistogramData
    label::String
    gene::String
    interval::Float64
    trace::traceType
    nRNA::Int
    histRNA::hType
end
