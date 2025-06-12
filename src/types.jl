"""
    StateDistTuple

Type alias for `Tuple{Int64, Vector{Distribution{Univariate,Continuous}}}`.
Represents a tuple containing a state index and a vector of continuous univariate probability distributions.
"""
const StateDistTuple = Tuple{Int64, Vector{Distribution{Univariate,Continuous}}} 