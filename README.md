# StochasticGene.jl

Julia package to simulate and fit stochastic models of gene transcription to experimental data. The data can range from distributions of mRNA counts per cell (either through single molecule FISH (smFISH) or single cell RNA sequence (scRNA) data to dwell time distributions of single pre-RNA molecules in the act of transcription imaged in live cells. The models are continuous Markov systems with an arbitrary number of G (gene) states and R (pre-RNA) steps with stochastic transitions between the states. The gene always occupies one of the G states and there are reversible transitions between G states.  One of the G states is an active state where transcription can be initiated and the first R step becomes occupied. An irreversible forward transition can then occur to the next R step if that step is unoccupied simulating elongation. An mRNA molecule is ejected from the final (termination) R step where it then decays stochastically. The model can account for multiple alleles of the gene in the same cell and coupling between alleles. Each R step is considered visible when occupied; the first R step represents the time the inserted reporter is first observable. In the original model in Rodriguez et al. Cell (2018), the reporter was in the exon and thus was carried out to the last step and ejected. In Wan et al. (submitted), the reporter is inserted into an intron and thus can be spliced out before the polymerase reaches the final R step. Models are allowed to have no R steps (i.e. classic telegraph models but with arbitrary numbers of G states) where an mRNA molecule can be stochastically produced when the gene is occupies the active G state.  When fitting the model to single cell RNA (scRNA) sequence data, the predicted mRNA counts are further subjected to a stochastic process to account for the reduction in the number of mRNA captured and measured.

The repository is organized into files that contain functions for specifying the models and preparing the data, computing the model predicted distributions for the dwell time distributions from live cell imaging measurements and intracellular mRNA distributions (either smFISH or scRNA), and applying a Bayesian Metropolis-Hastings markov chain monte carlo (MCMC) algorithm to fit the models to the data and compute posterior distributions of the parameters.

The functions to specify and prepare the data and models for three specific data sets have been provided. The models and data are organized into Julia structs. The MCMC alogorithm will call for function methods specific to the models and data for the prior and proposal distributions and the loglikelihood of the model prediction.

### Example:
To fit the four rates of a simple 2 G state model (i.e. telegraph model, G0<=>G1 -> mRNA -> 0) to smFISH mRNA data for the gene MYC,
you would launch Julia from the directory where the file StochasticGene.jl is located and type:

```
include("StochasticGene.jl")
data,model,options = StochasticGene.steadystate_rna(datafolder,name,gene,nsets,r,decayprior,G,nalleles,fittedparam,cv,maxtime,nsamples,temp)
fit,stats,waic=StochasticGene.run_mh(data,model,options)
```

where
`data = data structure`,
`model = model structure`,
`options = MCMC run parameters`,
`datafolder = "datafolder/"  (folder where data is stored)`,
`name = "steady state smFISH" (text to describe the run)`,
`gene = "MYC" (gene name)`,
`nsets = 1 (number )`,
`G = 2 (Int)`,
`nalleles = 2 (Int)`,
`fittedparam = [1;2;3;4] (vector of indices for which rate parameters to be fit)`, 
`cv = 0.05 (coefficient of variation in proposal distribution or a covariance matrix)`,
`maxruntime = 60. (in seconds)`,
`nsamples = 10000 (number of MCMC samples)`,
`temp = 100.  (MCMC temperature)`,
`r = vector containing intial parameter guess`,
`fit = structure containing MCMC ouputs`,
`stats = structure containing parameter statistics from MCMC samples`,
`waic = tuple containing Watanabe-Akaike Information Criterion (WAIC) and SE`

The function steadystate_rna prepares the necessary structures and run_mh runs the MCMC for a total of nsamples steps with a maximum clocktime of maxtime.
The data in `datafolder/` needs to be in the form of a histogram (number of cells for a given mRNA count).  These functions can also be put into a script (e.g. runscript.jl) and run from the Unix command line with `Julia runscript.jl`.
