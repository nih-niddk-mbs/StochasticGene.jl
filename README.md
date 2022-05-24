# StochasticGene.jl

Julia package to simulate and fit stochastic models of gene transcription to experimental data. The data can range from distributions of mRNA counts per cell (either through single molecule FISH (smFISH) or single cell RNA sequence (scRNA) data) to dwell time distributions of single pre-RNA molecules in the act of transcription imaged in live cells. The models are continuous Markov systems with an arbitrary number of G (gene) states and R (pre-RNA) steps with stochastic transitions between the states. The gene always occupies one of the G states and there are reversible transitions between G states.  One of the G states is an active state where transcription can be initiated and the first R step becomes occupied. An irreversible forward transition can then occur to the next R step if that step is unoccupied simulating elongation. An mRNA molecule is ejected from the final (termination) R step where it then decays stochastically. The model can account for multiple alleles of the gene in the same cell and coupling between alleles. Each R step is considered visible when occupied; the first R step represents the time the inserted reporter is first observable. In the original model in Rodriguez et al. Cell (2018), the reporter was in the exon and thus was carried out to the last step and ejected. In Wan et al. Cell (2021), the reporter is inserted into an intron and thus can be spliced out before the polymerase reaches the final R step. Models are allowed to have no R steps (i.e. classic telegraph models but with arbitrary numbers of G states) where an mRNA molecule can be stochastically produced when the gene occupies the active G state.  When fitting the model to single cell RNA (scRNA) sequence data, the predicted mRNA counts are further subjected to a stochastic process to account for the reduction in the number of mRNA captured and measured (i.e. technical noise).

The package has functions to specify the models, prepare the data, compute the predicted dwell time distributions from live cell imaging measurements and intracellular mRNA distributions (either smFISH or scRNA), and apply a Metropolis-Hastings markov chain monte carlo (MCMC) algorithm to fit the parameters of the models to the data and compute posterior distributions.

StochasticGene is designed to run on large data sets on a multiprocessor machine such as NIH Biowulf. Thus there are functions that generate swarm files to be submitted and process and analyze the results.



### Using StochasticGene


### Installation

StochasticGene is a registered Julia package.  To install on Biowulf, log on and at the prompt type:

```
[username@biowulf ~]$ sinteractive --mem=64G
```
This generates an interactive session

```
[username@biowulf ~]$ module load julialang
```
Loads Julia for use

```
[username@biowulf ~]$ julia - t 1
```
Starts Julia (with a single thread):

To install StochasticGene run the following in the Julia REPL:

```
julia> ] add StochasticGene
```

You can check if all tests pass by running (currently not working)

```
julia> ] test StochasticGene
```

Command "]" brings you into the Julia Package environment, "Ctrl C" gets out

StochasticGene will be updated periodically, to update on Julia type

```
julia> ] update StochasticGene
```

StochasticGene requires a specific directory structure where data are stored and results are saved.  At the top is the `root` folder (e.g. "scRNA" or "RNAfits") with subfolders `data` and `results`. Inside `data` are two more folders  `alleles` and `halflives`,  containing allele numbers and half lives, respectively.  The command `rna_setup` will create the folder structure. New allele numbers and halflives for new cells can be added directly to the folders.  The files should be csv format and have the form `[cell name]_alleles.csv` or `[cell name]_halflife.csv`.

```
julia> using StochasticGene

julia> rna_setup("scRNA")
```

or any other name you choose for the root directory.

### Example Use on Biowulf:

Fit the scRNA histogram in all the genes in folder called "data/HCT_testdata" (which should exist if you ran `setup`) on NIH Biowulf by running a swarmfile.


First move into the root directory you created and launch julia:

```
[username@biowulf ~]$ cd scRNA

[username@biowulf ~]$ julia - t 1

```

Create swarm files using the command in the JULIA repl:
```
julia> using StochasticGene

julia> makeswarm(["CENPL","MYC"],cell="HCT116",maxtime = 600.,nchains = 8,fish = false,cycle=true,G=2,conds = "MOCK",resultfolder ="HCT_scRNA",datafolder = "HCT116_testdata/",root = ".")
```

The genes are listed as a vector of strings. You only need to type `using StochasticGene` once per session.
This will create two files.  The first will end in `.swarm` and the second will end in `.jl`

To fit all the genes in the data folder use:
```
julia> using StochasticGene

julia> makeswarm(cell="HCT116",maxtime = 600.,nchains = 8,fish = false,cycle=true,G=2,conds = "MOCK",resultfolder ="HCT_scRNA",datafolder = "HCT116_testdata/",nsets=1,root = ".")
```

To exit julia type:

```
julia> exit()
```

To run the swarm file, type at the command line:

```
[username@biowulf ~]$ swarm -f fit_HCT116-scRNA-ss_MOCK_2.swarm --time 24:00:00 -t 8  -g 24 --merge-output --module julialang
```

This will submit a job into the Biowulf queue.  To check the status of your job type:

```
[username@biowulf ~]$ sjobs
```

When the job finishes, Biowulf will create new swarm files in your folder. The fit results will be saved in the folder `results/HCT_scRNA`.  There will be three result files for each gene and model.  The file names will have the form

`[filetype]_[label]_[condition]_[gene name]_[modeltype written as consecutive numbers GRS]_[number of alleles].txt`

`_` (underscore) is a reserved character and should not be used in any of the file field such as `[label]`.

filetypes are:

`rates`, which contain 4 lines for maximum likelihood rates, mean rates, median rates

`measures`, which contain information about the quality of the fits

`param-stats`, which contain detailed statistics of the parameter posteriors (the MCMC samples are not saved)

In our example the files `rates_HCT116-scRNA-ss_MOCK_CENPL_2_2.txt`,`measures_HCT116-scRNA-ss_MOCK_CENPL_2_2.txt`,`param-stats_HCT116-scRNA-ss_MOCK_CENPL_2_2.txt` will be produced

The output convention is that underscore `_` is used to separate the 4 fields of a result file and thus should not be used in any of the fields.

A data frame of the results can be constructed in Julia using the write_dataframes(resultfolder,datafolder) function, e.g.

```
julia> using StochasticGene

julia> write_dataframes("results/HCT_scRNAtest","data/HCT116_testdata")
```
The result will be a set of csv files collating the result files of all the genes along with a "Winners" file that lists which model performed best for a given measure (default is AIC but can be changed with a named argument, see API below) and a Summary file condensing the information of the other files.

The Summary file can be supplemented with more information via the write_augmented(summaryfile,resultfolder,datafolder) function, e.g.

```
julia> write_augmented("Summary_HCT116-scRNA-ss_MOCK_2.csv","results/HCT_scRNAtest","data/HCT116_testdata")


```

### Example Use on Unix
If not running on Biowulf, the same swarm files can be used, although they will not be run in parallel.

In Bash, type:

```
Bash> chmod 744 fit_scRNA-ss-MOCK_2.swarm

Bash> bash fit_scRNA-ss-MOCK_2.swarm &

```

This will execute each gene in the swarm file sequentially. To run several genes in parallel, you can break the run up into multiple swarm files and execute each swarm file separately.  You can also trade off time with the number of chains, so if you have a large processor machine run each gene on many processors, e.g. 64, for a short amount of time, e.g. 15 min.

### API:

```
makeswarm(;G::Int=2,cell="HCT116",swarmfile::String="fit",label="label",inlabel=label,timestamp="",nsets=1,datafolder::String="HCT116_testdata",fish= false,cycle=true,thresholdlow::Float64=0.,thresholdhigh::Float64=1e8,conds::String="MOCK",resultfolder::String= "fit_result",infolder=resultfolder,batchsize=1000,maxtime = 60.,nchains::Int=2,nthreads::Int=1,transient::Bool=false,fittedparam=collect(1:2*G-1),fixedeffects=(),juliafile::String="fitscript",root=".",samplesteps::Int=100000,warmupsteps=20000,annealsteps=0,temp=1.,tempanneal=100.,cv = 0.02,yieldprior=0.05)

makeswarm(genes::Vector;G::Int=2,cell="HCT116",swarmfile::String="fit",label="label",inlabel=label,timestamp="",nsets=1,datafolder::String="HCT116_testdata",fish=false,cycle=true,conds::String="MOCK",resultfolder::String="fit_result",infolder=resultfolder,batchsize=1000,maxtime=60.,nchains::Int=2,nthreads::Int=1,transient::Bool=false,fittedparam=collect(1:2*G-1),fixedeffects=(),juliafile::String="fitscript",root=".",samplesteps::Int=100000,warmupsteps=20000,annealsteps=0,temp=1.,tempanneal=100.,cv=0.02,yieldprior=0.05)

Arguments
    - `G`: number of gene states
    - `cell': cell type for halflives and allele numbers
    - `infolder`: name of folder for initial parameters
    - `swarmfile`: name of swarmfile to be executed by swarm
    - `label`: label of output files produced
    - `inlabel`: label of files used for initial conditions
    - `nsets`: number of histograms to be fit (e.g. one for wild type and one for perturbation)
    - `datafolder`: folder holding histograms, if two folders use `-` (hyphen) to separate, e.g.  "data\folder1-data\folder2"
    - `fish`: Data file type, set to true for FISH and false to scRNA (FISH type assumes different folder structure)
    - `thresholdlow`: lower threshold for halflife for genes to be fit
    - `threhsoldhigh`: upper threshold
    - `conds`: string describing conditions to be fit with `-` to separate if two conditions, e.g. "WT-AUXIN"
    - `result`: folder for results
    - `batchsize`: number of jobs per swarmfile, default = 1000
    - `maxtime`: maximum wall time for run, default = 2 hrs
    - `nchains`: number of MCMC chains = number of processors called by Julia, default = 2
    - 'nthreads`: number of Julia threads per processesor, default = 1
    - `transient::Bool`: true means fit transient model (T0, T30, T120)
    - `fittedparam`: vector of rate indices to be fit, e.g. [1,2,3,5,6,7]
    - `fixedeffects`: tuple of vectors of rates that are fixed between control and treatment where first index is fit and others are fixed to first, e.g. ([3,8],) means  index 8 is fixed to index 3
         (each vector in tuple is a fixed rate set)
    - `juliafile`: name of file to be called by julia in swarmfile
    - `root`: name of root directory for project, e.g. "scRNA\"
    - `samplesteps`: number of MCMC sampling steps
    - `warmupsteps`: number of MCMC warmup steps to find proposal distribution covariance
    - `annealsteps`: number of annealing steps (during annealing temperature is dropped from tempanneal to temp)
    - `temp`: MCMC temperature
    - `tempanneal`: annealing temperature
    - `cv`: coefficient of variation (mean/std) of proposal distribution
    - `yieldprior`: prior for yield, default = .05, (set to 1 for FISH if using scRNA data format for FISH data)
    - `genes`: array of genes to be fit


returns swarmfile that calls a julia file that is executed on biowulf

```

```
  write_dataframes(resultfolder::String,datafolder::String;measure::Symbol=:AIC,assemble::Bool=true)

  collates run results into a csv file

Arguments
- `resultfolder`: name of folder with result files
- `datafolder`: name of folder where data is stored
- `measure`: measure used to assess winner
- `assemble`: if true then assemble results into summary files
```

```
write_winners(resultfolder,measure)

Write best performing model for measure (e.g. AIC, WAIC, Deviance)

```

```
write_augmented(summaryfile::String,resultfolder,datafolder)

Augment summary file with burst size (for G > 1), model predicted moments, and fit measures


```
