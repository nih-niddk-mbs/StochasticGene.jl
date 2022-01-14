# StochasticGene.jl

Julia package to simulate and fit stochastic models of gene transcription to experimental data. The data can range from distributions of mRNA counts per cell (either through single molecule FISH (smFISH) or single cell RNA sequence (scRNA) data to dwell time distributions of single pre-RNA molecules in the act of transcription imaged in live cells. The models are continuous Markov systems with an arbitrary number of G (gene) states and R (pre-RNA) steps with stochastic transitions between the states. The gene always occupies one of the G states and there are reversible transitions between G states.  One of the G states is an active state where transcription can be initiated and the first R step becomes occupied. An irreversible forward transition can then occur to the next R step if that step is unoccupied simulating elongation. An mRNA molecule is ejected from the final (termination) R step where it then decays stochastically. The model can account for multiple alleles of the gene in the same cell and coupling between alleles. Each R step is considered visible when occupied; the first R step represents the time the inserted reporter is first observable. In the original model in Rodriguez et al. Cell (2018), the reporter was in the exon and thus was carried out to the last step and ejected. In Wan et al. (submitted), the reporter is inserted into an intron and thus can be spliced out before the polymerase reaches the final R step. Models are allowed to have no R steps (i.e. classic telegraph models but with arbitrary numbers of G states) where an mRNA molecule can be stochastically produced when the gene is occupies the active G state.  When fitting the model to single cell RNA (scRNA) sequence data, the predicted mRNA counts are further subjected to a stochastic process to account for the reduction in the number of mRNA captured and measured.

There are functions to specify the models, prepare the data, compute the model predicted distributions for the dwell time distributions from live cell imaging measurements and intracellular mRNA distributions (either smFISH or scRNA), and apply a Metropolis-Hastings markov chain monte carlo (MCMC) algorithm to fit the models to the data and compute posterior distributions of the parameters.

StochasticGene is designed to run on large data sets on a multiprocessor machine such as NIH Biowulf.



### Using StochasticGene


### Installation

StochasticGene is a registered Julia package.  To run Julia on Biowulf first log onto Biowulf and at the prompt type:

```
[username@biowulf ~]$ sinteractive --mem=64G
```

```
[username@biowulf ~]$ module load julialang
```

```
[username@biowulf ~]$ julia
```

To install StochasticGene run the following in the Julia REPL:

```
julia> ] add StochasticGene
```

You can check if all tests pass by running

```
julia> ] test StochasticGene
```

Command "]" brings you into the Julia Package environment, "Ctrl C" gets out

StochasticGene requires a specific directory structure where data are stored and results are saved.  At the top is the `root` folder (e.g. "scRNA" or "RNAfits") with subfolders `data` and `results`. Inside `data` are two more folders containing allele numbers and halflives.  The command `rna_setup` will create the folder structure.

```
julia> rna_setup("scRNA")
```

or any other name you choose for the root directory.

### Example Use:

Fit the scRNA histogram in all the genes in folder called "data/HCT_testdata" (which should exist if you ran `setup`) on NIH Biowulf by running a swarmfile

First create swarm files using the command in the JULIA repl:

```
makeswarm(maxtime = 600.,fish = false,cycle=true,G=2,conds = "MOCK",fittedparam=[1,2,3],resultfolder ="HCT_scRNA",datafolder = "HCT116_testdata/",root)
```

This will create two files.  The first will end in .swarm and the second will end in .jl

To exit julia type:

```
julia> exit()
```

To run the swarm file, type at the command line:

`
[username@biowulf ~]$ swarm -f fit_scRNA-ss-MOCK_3.swarm --time 24:00:00 -t 8  -g 24 --merge-output --module julialang
`

The results will be saved in the folder "results/HCT_scRNA".  There will be three files for each gene and model.  The file names will have the form

`[filetype]_[label]_[gene name]_[modeltype written as consecutive numbers GRS]_[number of alleles].txt`

filetypes are:

`rates`, which contain 4 lines for maximum likelihood rates, mean rates, median rates

`measures`, which contain information about the quality of the fits

`param_stats`, which contain detailed statistics of the parameter posteriors (the MCMC samples are not saved)

In our example the files `rates_FISH-ss-MOCK_CENPL_2_2.txt`,`measures_FISH-ss-MOCK_CENPL_2_2.txt`,`param-stats_FISH-ss-MOCK_CENPL_2_2.txt` will be produced

The output convention is that underscore `_` is used to separate the 4 attributes of a run and thus should not be used elsewhere.

A data frame of the results can be constructed in Julia using the commands

```
julia> using StochasticGene
```

```
julia> write_dataframe(rootfolder,outputfile,"HCT_FISHtest",["1","2"],"MOCK","data/HCT116_testdata",true)
```

###API
```
makeswarm(;G::Int=2,cell="HCT116",swarmfile::String="fit",label="label",inlabel=label,nsets=1,datafolder::String="data/HCT116_testdata",fish= false,cycle=true,thresholdlow::Float64=0.,thresholdhigh::Float64=1e8,conds::String="DMSO",resultfolder::String= "fit_result",infolder=resultfolder,batchsize=1000,maxtime = 60.,nchains::Int = 2,transient::Bool=false,fittedparam=[1],fixedeffects=(),juliafile::String="fitscript",root="../",samplesteps::Int=100000,warmupsteps=20000,annealsteps=0,temp=1.,tempanneal=100.,modulepath = "/Users/carsonc/github/StochasticGene/src/StochasticGene.jl",cv = 0.02)
```

```
makeswarm(genes::Vector;G::Int=2,cell="HCT116",swarmfile::String="fit",label="label",inlabel=label,nsets=1,datafolder::String="data/HCT116_testdata",fish=false,cycle=true,conds::String="DMSO",resultfolder::String="fit_result",infolder=resultfolder,batchsize=1000,maxtime=60.,nchains::Int=1,transient::Bool=false,fittedparam=[1],fixedeffects=(),juliafile::String="fitscript",root="../",samplesteps::Int=100000,warmupsteps=20000,annealsteps=0,temp=1.,tempanneal=100.,cv=0.02)
```

create swarm files to be run on NIH Biowulf

Arguments
- `G`: number of gene states
- `cell': cell type used for half life and allele number
- `infolder`: name of folder for initial parameters
- `swarmfile`: name of swarm file to be executed by swarm
- `label`: label of output files produced (e.g. "scRNA-ss-CONTROL")
- `inlabel`: label of files used for initial conditions
- `nsets`: number of histograms to be fit simultaneously (e.g. one for wild type and one for perturbation)
- `datafolder`: folder holding histograms, if two folders use `-` (hyphen) to separate, e.g.  "data\folder1-data\folder2"
- `thresholdlow`: lower threshold for half life for genes to be fit
- `threhsoldhigh`: upper threshold
- `conds`: string describing conditions to be fit with `-` to separate if two conditions, e.g. "WT-AUXIN"
- `result`: folder for results
- `batchsize`: maximum number of jobs per swarmfile, default = 1000
- `maxtime`: maximum wall time for run, default = 2 hrs
- `nchains`: number of MCMC chains, default = 2
- `transient::Bool`: true means fit transient model, e.g. at times 0, 30, and 120 minutes, using a time dependent model
- `fittedparam`: vector of rate indices to be fit, e.g. [1,2,3] means rate1, rate2 and rate3 are fit and all other rates are fixed
- `fixedeffects`: tuple of vectors of rates that are fixed between control and treatment where first index is fit and others are fixed to first, e.g. ([3,8],) means  index 8 is fixed to index 3 (each vector in tuple is a fixed rate set) only used for nsets > 1
- `juliafile`: name of file to be called by julia in swarmfile
- `root`: name of root directory for project, e.g. "scRNA\"
- `samplesteps`: number of MCMC sampling steps
- `warmupsteps`: number of MCMC warmup steps to find proposal distribution covariance
- `annealsteps`: number of annealing steps (during annealing temperature is dropped from tempanneal to temp)
- `temp`: MCMC temperature
- `tempanneal`: annealing temperature
- `cv`: coefficient of variation (mean/std) of proposal distribution

- `genes`: array of genes to be fit

```
write_dataframe(root::String,outputfile::String,folder::String,models::Vector,cond::String,datafolder::String,fish::Bool)
```
collates fit data into a date frame and saves into a csv file

Arguments
- `root`: root folder
- `outputfile`: name of an ouput file
- `folder`: name of folder with result files
- `models`: vector of models (listed as strings), e.g. ["1","2"]
- `cond`: condition of experiment
- `datafolder`: name of folder where data is stored
- `fish`: true if data is a FISH histogram (i.e. no technical loss is accounted for)
