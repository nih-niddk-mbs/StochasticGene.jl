# StochasticGene.jl

Julia package to simulate and fit stochastic models of gene transcription to experimental data. The data include distributions of mRNA counts per cell (e.g. single molecule FISH (smFISH) or single cell RNA sequence (scRNA) data)), image intensity traces from live cell imaging, and dwell time distributions of reporters (e.g. MS2 or PP7) imaged in live cells. The models are continuous Markov systems with an arbitrary number of G (gene) states, R (pre-RNA) steps with stochastic transitions between the states, and S splice sites (usually same as R). The gene always occupies one of the G states and there are reversible transitions between G states that you specify.  One of the G states is an active state where transcription can be initiated and the first R step becomes occupied. An irreversible forward transition can then occur to the next R step if that step is unoccupied simulating elongation. An mRNA molecule is ejected from the final (termination) R step where it then decays stochastically. The model can account for multiple alleles of the gene in the same cell. The user can specify which R step is considered visible when occupied. For example, if the pre-RNA MS2 construct is inserted into an intron that is transcribed early then you would select the reporter insertstep to be one and all R steps are visible. But if the reporter is transcribed late then you would choose a later R step, like step 3, and only R steps after step 3 are considered visible. In the original model in Rodriguez et al. Cell (2018), the reporter was in the exon and thus was visible to the last step and then ejected. In Wan et al. Cell (2021), the reporter is inserted into an intron and thus can be spliced out before the polymerase reaches the final R step. Models are allowed to have no R steps (i.e. classic telegraph models but with arbitrary numbers of G states) where an mRNA molecule can be stochastically produced when the gene occupies the active G state.  When fitting the model to single cell RNA (scRNA) sequence data, the predicted mRNA counts are further subjected to a stochastic process to account for the reduction in the number of mRNA captured and measured (i.e. technical noise). You can also specify reporters for G states and even more than one simultaneous reporter.

The package has functions to specify the models, prepare the data, compute the predicted data, apply a Metropolis-Hastings markov chain monte carlo (MCMC) algorithm to fit the parameters of the models to the data and compute posterior distributions, and simulate the models.

StochasticGene can run on small data sets on a laptop or large data sets on a multiprocessor cluster such as NIH Biowulf. There are functions that generate swarm files to be submitted and process and analyze the results.



### Using StochasticGene

### Installation

The following assumes that Julia has already been installed. If not go to https://julialang.org/downloads/. Julia has already been installed on the NIH Biowulf system but StochasticGene must be installed by each individual user.

To install StochasticGene, open a terminal and type

```
Bash> julia

After Julia opens you will be in the interactive Julia REPL.

To install StochasticGene run the following:

```
julia> ] add StochasticGene
```

The command "]" brings you into the Julia Package environment, "Ctrl C" gets you out

After the package has installed, you can check if it works by running

```
julia> ] test StochasticGene
```
(this could take 10 or more minutes)

StochasticGene will be updated on Github periodically, to update to a new version type

```
julia> ] update StochasticGene
```


To install StochasticGene on Biowulf, log on and at the prompt type:

```
[username@biowulf ~]$ sinteractive --mem=64G
```
This generates an interactive session

```
[username@biowulf ~]$ module load julialang
```
This loads Julia for use.  These two steps are only necessary on Biowulf. 

```
[username@biowulf ~]$ julia - t 1
```
Starts Julia (with a single thread). Then continue as before



StochasticGene requires a specific directory structure where data are stored and results are saved.  At the top is the `root` folder (e.g. "scRNA" or "RNAfits") with subfolders `data` and `results`. Inside `data` are two more folders  `alleles` and `halflives`,  containing allele numbers and half lives, respectively.  The command `rna_setup()` will create the folder structure. New allele numbers and halflives for new cells can be added directly to the folders.  These files should be csv format and have the form `[cell name]_alleles.csv` or `[cell name]_halflife.csv`.  Halflives for the then genes in Wan et al. for HBEC cells are built into the system.

After installing StochasticGene, you can set up the folder structure by typing

```
julia> using StochasticGene

julia> rna_setup("scRNA")
```
or any other name you choose for the root directory.

### Fitting data with StochasticGene

To fit a model, you need to load data and choose a model. Data allowed are stationary histograms (e.g. smFISH or scRNA), intensity traces (e.g. trk files), nascent RNA fraction (e.g. from intronic FISH probes), and dwell time distributions. Different data types from the same experiment can also be fit simultaneously. For example, rna histograms from smFISH can be fit together with multiple traces or dwell time histograms from live cell recordings.

Models are distringuished by the number of G states, transitions between G states, number of R steps, the step where the reporter is inserted, and whether splicing occurs.  For intensity traces and dwell time distributions, the sojourn states (i.e. R steps or G states where the reporter is visible) called "on states" must also be specified.

Data can be fit using the `fit` function, which has named arguments (see API below) to determine the type of data to be fit, where it can be found, what model to be used, and what options to be set for the fitting algorithm. The default, run by typing ```fit()```, will fit the mock rna histogram data installed by rna_setup(root) with a simple two state telegraph model (2 G stsates, no R steps).  Thwo files are in the folder "root/data/HCT116_testdata", where root is specified by the user. The default root is the folder in which julia was launched but can be specified using the named argument `root`. The `fit` function returns six variables, resulting in

```
julia> fits, stats, measures, data, model, options = fit();
2023-10-24T10:22:10.782
Gene: MYC G: 2 Treatment:  MOCK
data: HCT116_testdata/
in: HCT116_test out: HCT116_test
maxtime: 60.0
./results/HCT116_test/rates_rna-HCT116_MOCK_MYC_2_2.txt does not exist
[0.01, 0.01, 0.1, 0.032430525886402085]
initial ll: 38980.680174070265
final max ll: 21755.280380669064
median ll: 21755.280645984076
Median fitted rates: [0.020960969159742875, 0.142016620787, 1.1015638965800167]
ML rates: [0.020979229436880756, 0.1423818501792811, 1.10350961139495]
Acceptance: 325538/613219
Deviance: 0.007326880506251085
rhat: 1.0105375337987144
2023-10-24T10:24:11.295
```
The semicolon is not necessary but suppresses printing the returned variables.

The `fit` function prints out information about the time, the data, the model, and information about the results. A more detailed version is given in the returned variables and written to the folder "./results/HCT116_test". You can use the results of a previous run as an initial condition. In this case, there were no previous runs and thus the function used a default starting point, which is also the prior. In this particular run, only three rates were fitted, and the median of the posterior and the maximum likelihood parameters are printed out. `Acceptance` shows that in the allotted 60 seconds of real time, out of 613219 samples, 325538 were accepted by the Metropolis-Hastings MCMC algorithm. The `Deviance` is the difference in likelihood between a perfect fit and the given fit. Zero indicates a perfect fit and anything less than one is a good fit. `rhat` is a measure of how close to convergence the run achieved with 1 being ideal. In this particular run, only one chain was used. More chains can be specified and each run on it's own processor. To use more chains, you need to specify more processors by staring Julia with

```
Bash> julia -p 4

julia> @everywhere using StochasticGene

julia> fits, stats, measures, data, model, options = fit(nchains=4);
2023-10-24T10:55:16.232
Gene: MYC G: 2 Treatment:  MOCK
data: HCT116_testdata/
in: HCT116_test out: HCT116_test
maxtime: 60.0
[0.020960969159742875, 0.142016620787, 1.1015638965800167, 0.032430525886402085]
initial ll: 21755.280645984076
final max ll: 21755.280378544387
median ll: 21755.282721546962
Median fitted rates: [0.020936177089260818, 0.14140674439592604, 1.0983188350949695]
ML rates: [0.020970183926180535, 0.14229910450327884, 1.10310357519155]
Acceptance: 713078/1342174
Deviance: 0.0073268798846403485
rhat: 1.0018023969479748
2023-10-24T10:56:31.337

```


"rna", "rnaonoff", "rnadwelltime", "trace", "tracenascent", and "tracerna".

```


 The fit function.


a=fit(4,"rnaoffon",[],["FISH","dwelltime"],"DNAJC5","HBEC","",5/3,0.,"10-4","10-4","gt","gt",[1, 2, 3, 4],tuple(),transitions,2,1,0,1,"/Users/carsonc/Box/Larson/GeneTrap_analysis/")

a=fit(4,"rnatrace",[],["FISH","traces"],"DNAJC5","HBEC","",5/3,0.,"10-4","10-4","gt","gt",[1, 2, 3, 4, 6,7,8,9,10],tuple(),transitions,2,1,0,1,"/Users/carsonc/Box/Larson/GeneTrap_analysis/")

fit(4,"rnadwelltime",[],["FISH","dwelltime"],"DNAJC5","HBEC","",5/3,0.,"10-4","10-4","gt","gt",[1, 2, 3, 4, 6,7,8,9,10],tuple(),transitions,2,1,0,1,"/Users/carsonc/Box/Larson/GeneTrap_analysis/")


### Example Use on Unix
If not running on Biowulf, the same swarm files can be used, although they will not be run in parallel.

In Bash, type:

```
Bash> chmod 744 fit_scRNA-ss-MOCK_2.swarm

Bash> bash fit_scRNA-ss-MOCK_2.swarm &

```

This will execute each gene in the swarm file sequentially. To run several genes in parallel, you can break the run up into multiple swarm files and execute each swarm file separately.  You can also trade off time with the number of chains, so if you have a large processor machine run each gene on many processors, e.g. 64, for a short amount of time, e.g. 15 min.

### Example Use on Biowulf to fit scRNA data:


Fit the scRNA histogram in all the genes in folder called "data/HCT_testdata" (which should exist if you ran `setup`) on NIH Biowulf by running a swarmfile.


First move into the root directory you created and launch julia:

```
[username@biowulf ~]$ cd scRNA

[username@biowulf ~]$ julia - t 1

```

Create swarm files using the command in the JULIA repl:
```
julia> using StochasticGene

julia> makeswarm(["CENPL","MYC"],cell="HCT116",maxtime = 600.,nchains = 8,datatype = "rna",G=2,transitions = ([1,2],[2,1]),datacond = "MOCK",resultfolder ="HCT_scRNA",datapath = "HCT116_testdata/",root = ".")
```

The genes are listed as a vector of strings. You only need to type `using StochasticGene` once per session.
This will create two files.  The first will end in `.swarm` and the second will end in `.jl`

To fit all the genes in the data folder use:
```
julia> using StochasticGene

julia> makeswarm(cell="HCT116",maxtime = 600.,nchains = 8,datatype = "rna",G=2,transitions = ([1,2],[2,1]), datacond = "MOCK",resultfolder ="HCT_scRNA",datapath = "HCT116_testdata/",nsets=1,root = ".")
```

(for G = 3, use transitions = ([1,2],[2,1],[2,3],[3,2]), for 3 state Kinetic Proofreading model use transitions = ([1,2],[2,1],[2,3],[3,1]))

To exit julia type:

```
julia> exit()
```

To run the swarm file, type at the command line:

```
[username@biowulf ~]$ swarm -f fit_HCT116-scRNA-ss_MOCK_2.swarm --time 12:00:00 -t 8  -g 24 --merge-output --module julialang
```
(choose a time longer than maxtime (remember to convert seconds to hours))

This will submit a job into the Biowulf queue.  To check the status of your job type:

```
[username@biowulf ~]$ sjobs
```

When the job finishes, Biowulf will create new swarm files in your folder. The fit results will be saved in the folder `results/HCT_scRNA`.  There will be three result files for each gene and model.  The file names will have the form

`[filetype]_[label]_[condition]_[gene name]_[modeltype written as consecutive numbers GRS]_[number of alleles].txt`

`_` (underscore) is a reserved character and should not be used in any of the file field such as `[label]`.

filetypes are:

`rates`, 4 lines: maximum likelihood rates, mean rates, median rates

`measures`, information about the quality of the fits

`param-stats`, detailed statistics of the parameter posteriors (the MCMC samples are not saved)

`burst`, statistics on burst frequency, 7 lines: mean, sd, median, mad, quantiles at: .025,.5,.97.5

`optimized`, 3 lines: LBFGS optimized rates, negative log likelihood, convergence

In our example the files `rates_HCT116-scRNA-ss_MOCK_CENPL_2_2.txt`,`measures_HCT116-scRNA-ss_MOCK_CENPL_2_2.txt`,`param-stats_HCT116-scRNA-ss_MOCK_CENPL_2_2.txt`, `burst_HCT116-scRNA-ss_MOCK_CENPL_2_2.txt`, `optimized_HCT116-scRNA-ss_MOCK_CENPL_2_2.txt` will be produced

The output convention is that underscore `_` is used to separate the 4 fields of a result file and thus should not be used in any of the fields.

A data frame of the results can be constructed in Julia using the write_dataframes(resultfolder,datafolder) function, e.g.

```
julia> using StochasticGene

julia> write_dataframes_only("results/HCT_scRNAtest","data/HCT116_testdata")
```
The result will be a set of csv files collating the result files of all the genes along with a "Winners" file that lists which model performed best for a given measure (default is AIC but can be changed with a named argument, see API below) and a Summary file condensing the information of the other files.

The Summary file can be supplemented with more information via the write_augmented(summaryfile,resultfolder) function, e.g.

```
julia> write_augmented("results/HCT_scRNAtest/Summary_HCT116-scRNA-ss_MOCK_2.csv","results/HCT_scRNAtest")


```

### Simulations
Simulate any GRSM model using function simulator, which can produce steady state mRNA histograms, simulated intensity traces, and ON and OFF live cell histograms as selected.

### Units
StochasticGene assumes all rates have units of inverse minutes and the half lives in the `halflives` file are in hours. When computing or fitting stationary mRNA distributions, the rate units are relative. Scaling all the rates by a constant will not affect the results. In these cases, it is sometimes convenient to scale all the rates by the mRNA decay time, which is the last entry of the rate array. The rate units matter when considering or evaluating traces and histograms of ON and OFF times. The code assumes that these dwell time histogram have units of minutes (i.e. the reciprocal of the rate units). 

### API:

```
 makeswarm(;G::Int=2,R::Int=0,S::Int=0,transitions=([1,2],[2,1]),cell="HCT116",swarmfile::String="fit",label="label",inlabel="label",timestamp="",nsets=1,datafolder::String="HCT116_testdata",datatype="scRNA",thresholdlow::Float64=0.,thresholdhigh::Float64=1e8,conds::String="MOCK",resultfolder::String= "fit_result",infolder=resultfolder,batchsize=1000,maxtime = 60.,nchains::Int=2,nthreads::Int=1,transient::Bool=false,fittedparam=collect(1:2*G-1),fixedeffects=(),juliafile::String="fitscript",root=".",samplesteps::Int=40000,warmupsteps=20000,annealsteps=0,temp=1.,tempanneal=100.,cv = 0.02,priorcv= 10.,decayrate=-1.,burst=true,nalleles=2,optimize=true,type="",rtype="median")

    function makeswarm(genes::Vector;G::Int=2,R::Int=0,S::Int=0,transitions=([1,2],[2,1]),cell="HCT116",swarmfile::String="fit",label="label",inlabel="label",timestamp="",nsets=1,datafolder::String="HCT116_testdata",datatype="scRNA",conds::String="MOCK",resultfolder::String="fit_result",infolder=resultfolder,batchsize=1000,maxtime=60.,nchains::Int=2,nthreads::Int=1,transient::Bool=false,fittedparam=collect(1:2*G-1),fixedeffects=(),juliafile::String="fitscript",root=".",samplesteps::Int=40000,warmupsteps=20000,annealsteps=0,temp=1.,tempanneal=100.,cv=0.02,priorcv=10.,decayrate=-1.,burst=true,nalleles=2,optimize=true,type="",rtype="median")

    Arguments
    - `G`: number of gene states
    - `R`: number of pre-RNA steps (set to 0 for classic telegraph models)
    - `S`: number of splice sites (set to 0 for classic telegraph models and R for GRS models)
    - `transitions`: tuple of vectors that specify state transitions for G states, e.g. ([1,2],[2,1]) for classic 2 state telegraph model and ([1,2],[2,1],[2,3],[3,1]) for 3 state kinetic proof reading model
    - `cell': cell type for halflives and allele numbers
    - `swarmfile`: name of swarmfile to be executed by swarm
    - `label`: label of output files produced
    - `inlabel`: label of files used for initial conditions
    - `timestamp`: label for time of sample (e.g. T120)
    - `nsets`: number of histograms to be fit (e.g. one for wild type and one for perturbation)
    - `datafolder`: folder holding histograms, if two folders use `-` (hyphen) to separate, e.g.  "data\folder1-data\folder2"
    - `datatype`: String that desecribes data file type, e.g. "scRNA", "fish", "genetrap"
    - `thresholdlow`: lower threshold for halflife for genes to be fit
    - `threhsoldhigh`: upper threshold
    - `conds`: string describing data treatment condition, e.g. "WT", "DMSO", use `-` for two conditions, e.g. "WT-AUXIN"
    - `resultfolder`: folder for results
    - `infolder`: folder for initial parameters
    - `batchsize`: number of jobs per swarmfile, default = 1000
    - `maxtime`: maximum wall time for run, default = 2 hrs
    - `nchains`: number of MCMC chains = number of processors called by Julia, default = 2
    - 'nthreads`: number of Julia threads per processesor, default = 1
    - `transient::Bool`: true means fit a time dependent transient model (T0, T30, T120)
    - `fittedparam`: vector of rate indices to be fit, e.g. [1,2,3,5,6,7]
    - `fixedeffects`: tuple of vectors of rates that are fixed between control and treatment where first index is fit and others are fixed to first, e.g. ([3,8],) means  index 8 is fixed to index 3
         (each vector in tuple is a fixed rate set)
    - `juliafile`: name of file to be called by julia in swarmfile
    - `fitscript`: name of fit file called by swarm file
    - `root`: name of root directory for project, e.g. "scRNA\"
    - `samplesteps`: number of MCMC sampling steps
    - `warmupsteps`: number of MCMC warmup steps to find proposal distribution covariance
    - `annealsteps`: number of annealing steps (during annealing temperature is dropped from tempanneal to temp)
    Arguments below usually do not need to be set
    - `temp`: MCMC temperature
    - `tempanneal`: annealing temperature
    - `cv`: coefficient of variation (mean/std) of proposal distribution, if cv <= 0. then cv from previous run will be used
    - 'priorcv`: coefficient of variation for the rate prior distributions, default is 10.
    - `decayrate`: decay rate of mRNA, if set to -1, value in decayrate file will be used
    - `burst`: if true then compute burst frequency
    - `nalleles`: number of alleles, it will be overridden by allele file if it exists
    - `optimize`: use optimizer to compute maximum likelihood value
    - `type`: switch used for GRS models, choices include "", "offeject"
    - `rtype`: which rate to use for initial condition, choices are "ml", "mean", "median", or "last"
    - `writesamples`: write out MH samples if true, default is false


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
write_augmented(summaryfile::String,resultfolder::String)

Augment summary file with burst size (for G > 1), model predicted moments, and fit measures


```

```
	simulator(r::Vector{Float64},transitions,G::Int,R::Int,S::Int,nhist::Int,nalleles::Int;onstates::Vector{Int}=[G],range::Vector{Float64}=Float64[],total::Int=10000000,tol::Float64=1e-6,verbose::Bool=false)

	Simulate any GRSM model. Returns steady state mRNA histogram and if range not a null vector will return ON and OFF time histograms.
    If trace is set to true, it returns a nascent mRNA trace

	Arguments
	  - `r`: vector of rates
	  - `transitions`: tuple of vectors that specify state transitions for G states, e.g. ([1,2],[2,1]) for classic 2 state telegraph model and 
    ([1,2],[2,1],[2,3],[3,1]) for 3 state kinetic proof reading model
	  - `G`: number of gene states
    - `R`: number of pre-RNA steps (set to 0 for classic telegraph models)
    - `S`: number of splice sites (set to 0 for G (classic telegraph) and GR models and R for GRS models)
	  - `nhist::Int`: Size of mRNA histogram
	  - `nalleles`: Number of alleles

	Named arguments
    - `onstates::Vector`: a vector of ON G states
	  - `range::Vector{Float64}=Float64[]`: vector of time bins for ON and OFF histograms
	  - `totalsteps::Int=10000000`: maximum number of simulation steps
	  - `tol::Float64=1e-6`: convergence error tolerance for mRNA histogram (not used when outputting traces)
    - `traceinterval`: Interval in minutes between frames for intensity traces.  If 0, traces are not made. May need to try different totalsteps to get
    desired length of trace
    - `verbose::Bool=false`: flag for printing state information


  ### Examples

  julia> trace = simulator([.1,.02,.1,.05,.01,.01],([1,2],[2,1],[2,3],[3,1]),3,0,0,100,1,onstates=[2,3],traceinterval=100.,totalsteps = 1000)

  julia> hoff,hon,mhist = simulator([.1,.02,.1,.05,.01,.01],([1,2],[2,1],[2,3],[3,1]),3,0,0,20,1,onstates=[2,3],range=collect(1.:100.))


```

```
new_FISH(newroot::String,oldroot::String,rep::String)

Create new folder for FISH data with only one replicate

Arguments
  - `newroot`: new data folder e.g. HCT_T120_FISH_rep1
  - `oldroot`: old data folder e.g. HCT_T120_FISH
  - `rep`: name of replicate folder e.g. rep1

### Example

julia> new_FISH("data/HCT_T120_8genes_FISH_rep1","data/HCT_T120_8genes_FISH","rep1")


```
