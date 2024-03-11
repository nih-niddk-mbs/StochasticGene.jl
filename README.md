# StochasticGene.jl

Julia package to simulate and fit stochastic models of gene transcription to experimental data. The data include distributions of mRNA counts per cell (e.g. single molecule FISH (smFISH) or single cell RNA sequence (scRNA) data)), image intensity traces from live cell imaging, and dwell time distributions of reporters (e.g. MS2 or PP7) imaged in live cells. The transcription models considered are stochastic (continuous Markov systems) with an arbitrary number of G (gene) states, R (pre-RNA) steps, and S splice sites (usually same as R). The gene always occupies one of the G states and there are random (user specified) transitions between G states.  One of the G states is an active state where transcription can be initiated and the first R step becomes occupied. An irreversible forward transition can then occur to the next R step if that step is unoccupied mimicking elongation. An mRNA molecule is ejected from the final (termination) R step where it then decays stochastically. The model can account for multiple alleles of the gene in the same cell. The user can specify which R step is considered visible when occupied (i.e. contains a reporter). For example, if the pre-RNA MS2 construct is inserted into an intron that is transcribed early then you could make the reporter insertstep to be 1 and all R steps will be visible. But if the reporter is transcribed late then you could choose a later R step, like step 3, and only R steps after step 3 are considered visible. In the original model in Rodriguez et al. Cell (2018), the reporter was in an exon and was visible at all steps and then ejected. In Wan et al. Cell (2021), the reporter is inserted into an intron and thus can be spliced out before the polymerase reaches the final R step. Models are allowed to have no R steps (i.e. classic telegraph models but with arbitrary numbers of G states (rather thabn usual two)) where an mRNA molecule can be stochastically produced when the gene occupies the active G state. You can also specify reporters for G states (e.g. transcription factors) and more than one simultaneous reporter (e.g. reporters for introns and transcription factors).

The package has functions to specify the models, prepare the data, compute the predicted data, apply a Metropolis-Hastings markov chain monte carlo (MCMC) algorithm to fit the parameters of the models to the data (i.e. compute Bayesian posterior distributions), and explicitly simulate models.

StochasticGene can run on small data sets on a laptop or large data sets on a multiprocessor cluster such as NIH Biowulf. There are functions that generate swarm files to be submitted and process and analyze the results.


### Installing StochasticGene

The following assumes that Julia has already been installed. If not go to https://julialang.org/downloads/. Julia has already been installed on the NIH Biowulf system but StochasticGene must be installed by each individual user.

To install StochasticGene on a local computer, open a terminal and type

```
$ julia
```

After Julia opens you will be in the interactive Julia REPL.

Run the following:

```
julia> ] add StochasticGene
```

The command "]" brings you into the Julia Package environment, "Ctrl C" gets you out

After the package has installed, you can check if it works by running

```
julia> ] test StochasticGene
```
(this could take 10 or more minutes)

StochasticGene will be updated on Github periodically. To update to a new version type

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
This loads Julia for use. 

```
[username@biowulf ~]$ julia - t 1
```
Starts Julia (with a single thread). Then continue as before

Note: Sometimes Julia has problems crashing on Biowulf if an update is made to StochasticGene.  The best way to deal with this is to go to your
home directory and delete the julia directory via:
```
[username@biowulf ~]$ rm -r --force .julia
```
Then start julia and add StochasticGene again.

### Creating folder structure to run StochasticGene

StochasticGene requires a specific directory or folder structure where data are stored and results are saved.  At the top is the `root` folder (e.g. "scRNA" or "RNAfits") with subfolders `data` and `results`. Inside `data` are two more folders  `alleles` and `halflives`,  containing allele numbers and half lives, respectively.  The command 

```
julia> using StochasticGene

julia> rna_setup("scRNA")
```
will create the folder structure in the folder `scRNA` with the `data` and `results` subfolders, with some mock data in the `data` folder. Typing `rna_setup()` will assume the current folder is the root. Files for allele numbers and halflives for desired cell types can be added directly to the `data` folder.  These files should be csv format and have the form `[cell name]_alleles.csv` or `[cell name]_halflife.csv`.  Halflives for the then genes in Wan et al. for HBEC cells are built into the system.

### Fitting data with StochasticGene

To fit a model, you need to load data and choose a model. Data currently accepted are stationary histograms (e.g. smFISH or scRNA), intensity time traces (e.g. trk files), nascent RNA data (e.g. from intronic FISH probes), and dwell time distributions. Different data types from the same experiment can also be fit simultaneously. For example, RNA histograms from smFISH can be fit together with multiple intensity traces or dwell time histograms from live cell recordings.

Models are distringuished by the number of G states, the transition graph between G states (i.e. which G states can transitions to which other G states), number of R steps, the step where the reporter is inserted, and whether splicing occurs.  For intensity traces and dwell time distributions, the sojourn or "on" states (i.e. R steps or G states where the reporter is visible) must also be specified. Multiple sets of on states are allowed.

Data are fit using the `fit` function, which has named arguments (see API below) to determine the type of data to be fit, where it can be found, what model to be used, and what options to be set for the fitting algorithm. The default, run by typing `fit()`, will fit the mock rna histogram data installed by `rna_setup(root)` with a simple two state telegraph model (2 G stsates, no R steps).  Data is in the folder "root/data/HCT116_testdata", where root is specified by the user. The default root is the folder in which julia was launched but can be specified using the named argument `root`. The `fit` function returns six variables. 

To fit on a single processor, type

```
julia> using StochasticGene

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
The semicolon is not necessary but suppresses printing the returned variables. You only need to type `using StochasticGene` once per session.

The data is fit using the function `run_mh(data,model,options,nchains)` (which runs a Metropolis-Hastings MCMC algorithm), where `data`, `model`, and `options` are structures (Julia types) constructed by `fit` and returned. `nchains` is the number of MCMC chains, each running on a separate processor.  The three structures `fits`, `stats`, and `measures` (returned by `run_mh` and `fit`, give the fit results (Bayesian posteriors) and measures of the quality of the fit and are also saved to the folder "./results/HCT116_test", which is specified by the `resultfolder` argument.

During the run, the `fit` function prints out some of the information in fits, stats, and measures structures. `fit` will look for previous results in the folder specified by the argument `infolder` as an initial condition to start the MCMC run. In this case, there were no previous runs and thus the default was used, which is the priors of the parameters. In this run three rates were fit (set by argument `fittedparams`); the median of the posterior and the maximum likelihood parameters of the resulting fits are printed. The run was capped to a maximum of 60 seconds of real time (set by the argument `maxtime`) and `Acceptance` shows that out of 613219 MCMC samples, 325538 were accepted by the MCMC algorithm. The positive number `Deviance` is the difference in likelihood between a perfect fit and the resulting fit (for histograms only). `rhat` is a measure of MCMC convergence with 1 being ideal. 

The two state telegraph model has 4 transition rates, which are stored in a single vector. The order is 1) G state transitions (in the order as specified by the transition vector), 2) eject rate, 3) decay rate.  For general GRSM models, the order is 1) G state transitions 2) R transitions, 3) S transitions, 4) decay.  If there is no splicing then the S transitions are left out.  Not all the rates need to be fitted. In the above example, the decay rate is not fitted. This is specified by the fittedparams argument (see API). The posterior median, mean, standard deviations, and MADs will only be for the fitted params.

In this particular run, only one chain was used so the measures are not very informative. To use more chains, specify more processors with

```
$ julia -p 4

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

The `datatype` argument is a String that specifies the types of data to be fit. Currently, there are six choices: 1) "rna", which expects a single rna histogram file where the first column is the histogram; 2) "rnaonoff", which expects a rna histogram file and a three column dwelltime file with columns: bins, ON time distribution, and OFF time distribution; 3) "rnadwelltime", which fits an rna histogram together with multiple dwell time histograms specified by a vector of dwell time types, the choices being "ON","OFF" for R step reporters and "ONG", "OFFG" for G state reporters. Each dwelltime histogram is a two column file of the bins and dwell times and datapath is a vector of the paths to each; 4) "trace", which expects a folder of trace intensity (e.g. trk) files; 5) "tracenascent", the same with the nascent RNA input through the argument `nascent`; 6) "tracerna": trace files with RNA histogram. The data files are specified by the argument `datapath`. This can be a string pointing to a file or a folder containing the file.  If it points to a folder then `fit` will use the arguments `gene` and `datacond` to identify the correct file. For datatypes that require more than one data file, `datapath` is a vector of paths.  The path can include or not include the root folder.  

### Example fitting traces

Start Julia with multiple processors.

```
$ julia -p 4

```
Using `sinteractive` on Biowulf, you may need to also declare a single thread with
```
$ julia -p 4 -t 1

```
You can see how many processors have been called with

```

julia> nworkers()
4

julia @everywere using StochasticGene
```

You can create some simulated mock trace data with

```

julia> simulate_trace_data("data/testtraces/")

```

which generates 10 mock trk files in folder "data/testraces". See API below for the parameters used in the simulated data.

```

julia> fits, stats, measures, data, model, options = fit(datatype="trace",nchains=4,transitions=([1, 2], [2, 1], [2, 3], [3, 1]), G=3, R=2, S=2, insertstep=1, datapath="data/testtraces", gene="test", datacond="", decayrate=1.);
2023-10-28T15:39:16.955
Gene: test G R S insertstep: 3221
in: HCT116_test out: HCT116_test
maxtime: 60.0
./results/HCT116_test/rates_trace-HCT116__test_3221_2.txt does not exist
[0.01, 0.01, 0.01, 0.01, 0.1, 0.1, 0.1, 0.1, 0.1, 1.0, 100.0, 100.0, 100.0, 100.0, 0.9]
initial ll: 59437.54504442659
final max ll: 52085.58823772918
median ll: 55631.823149864926
Median fitted rates: [0.010164892707601266, 0.010415494521142342, 0.009888467865883247, 0.010565680996634689, 0.1051981476344431, 0.09544095249917071, 0.09473989790249088, 0.10045726361481444, 0.10418638776390472, 63.19680508738758, 56.794573780192, 119.72050615415071, 97.22874914820788, 0.9031904231186259]
ML rates: [0.008648315183072212, 0.01105818085890487, 0.01104728099474845, 0.011437834075767405, 0.09527550344438158, 0.11007539741898248, 0.08145019422061348, 0.09898368400856507, 0.1051451404740429, 37.71078406511133, 29.49680612047366, 171.40544543928038, 100.94919964017241, 0.9036456980377264]
Acceptance: 770/1525
rhat: 2.9216933196063533
2023-10-28T15:40:18.448

```
In this run, `rhat` is close to 3 indicating that the number of samples was insufficient to obtain a good sampling of the posterior distributions. either `maxtime` or `samplesteps` needs to be increased.


### Batch fitting on Biowulf using `swarm`.

The data can also be fit in batch mode using swarm files.  To do so, you create swarm files.

For example, you can fit all the scRNA histogram of all the genes in folder called "data/HCT_testdata" (which should exist if you ran `setup`) on NIH Biowulf by running a swarmfile.

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

### Example Use on Unix
If not running on Biowulf, the same swarm files can be used, although they will not be run in parallel.

e.g. in UNIX Bash, type:

```
Bash> chmod 744 fit_scRNA-ss-MOCK_2.swarm

Bash> bash fit_scRNA-ss-MOCK_2.swarm &

```

This will execute each gene in the swarm file sequentially. To run several genes in parallel, you can break the run up into multiple swarm files and execute each swarm file separately.  You can also trade off time with the number of chains, so if you have a large processor machine run each gene on many processors, e.g. 64, for a short amount of time, e.g. 15 min.

### Simulations
Simulate any GRSM model using function simulator, which can produce steady state mRNA histograms, simulated intensity traces, and ON and OFF live cell histograms as selected.

The following will simulate the steady state mRNA histogram, and ON and OFF dwelltime distributions, for an intron reporter inserted at the first R step, and for a transcription factor reporter visible in the G states 2 and 3.

```
julia> r=[0.038, 1.0, 0.23, 0.02, 0.25, 0.17, 0.02, 0.06, 0.02, 0.000231]
10-element Vector{Float64}:
 0.038
 1.0
 0.23
 0.02
 0.25
 0.17
 0.02
 0.06
 0.02
 0.000231

julia> transitions=([1, 2], [2, 1], [2, 3], [3, 1])
([1, 2], [2, 1], [2, 3], [3, 1])

julia> h=simulator(r,transitions,3,2,2,1,nhist=150,bins=[collect(5/3:5/3:200),collect(.1:.1:20)],onstates=[Int[],[2,3]],nalleles=2)
5-element Vector{Vector}:
 [102.03717039741856, 142.59069980711823, 88.32384512491423, 5.3810781196239645, 11.34932791013432, 20.045265169892332, 98.28644192386656, 58.86668672945706, 84.28930404506343, 38.33804408987862  …  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
 [729817, 620509, 530022, 456147, 397358, 349669, 315092, 284930, 261379, 241602  …  361, 369, 296, 300, 285, 266, 241, 216, 189, 190]
 [242326, 290704, 293370, 275618, 249325, 221256, 195222, 171730, 151167, 134806  …  18965, 18726, 18406, 18488, 18409, 18305, 17786, 17447, 17398, 17380]
 [2895640, 2557764, 2264680, 2001941, 1771559, 1567632, 1391091, 1228585, 1086723, 964156  …  7995, 7824, 7987, 7858, 8022, 7839, 7848, 7858, 7974, 7920]
 [116520, 116748, 115010, 115952, 114249, 114426, 114428, 113435, 113013, 112732  …  56077, 56462, 55967, 55863, 55710, 55549, 55092, 54988, 54984, 54528]

```

h is a vector containing 5 histograms. Although, not all data will include both ON and OFF time histograms, `simulator` automatically computes both.

### Units
StochasticGene assumes all rates have units of inverse minutes and the half lives in the `halflives` file are in hours. When computing or fitting stationary mRNA distributions, the rate units are relative. Scaling all the rates by a constant will not affect the results. In these cases, it is sometimes convenient to scale all the rates by the mRNA decay time, which is the last entry of the rate array. The rate units matter when considering or evaluating traces and histograms of ON and OFF times. The code assumes that these dwell time histogram have units of minutes (i.e. the reciprocal of the rate units). 

### API:

```
    fit(; <keyword arguments> )

Fit steady state or transient GM model to RNA data for a single gene, write the result (through function finalize), and return nothing.

#Arguments
- `nchains::Int=2`: number of MCMC chains = number of processors called by Julia, default = 2
- `datatype::String=""`: String that desecribes data type, choices are "rna", "rnaonoff", "rnadwelltime", "trace", "tracenascent", "tracerna"
- `dttype=String[]`
- `datapath=""`: path to data file or folder or array of files or folders
- `cell::String=""': cell type for halflives and allele numbers
- `datacond=""`: string or vector of strings describing data treatment condition, e.g. "WT", "DMSO" or ["DMSO","AUXIN"]
- `traceinfo=(1.0, 1., .65)`: 3 tuple of frame interval of intensity traces, transient (burn) time, and fraction of active traces
- `nascent=(1, 2)`: 2 tuple of number of spots, total number of locations (e.g. number of cells times number of alleles/cell)
- `traceinfo=tuple()`: tuple of trace information (frame interval, transient::Float64,onfraction::Float64)
- `infolder::String=""`: result folder used for initial parameters
- `resultfolder::String=test`: folder for results of MCMC run
- `label::String=""`: label of output files produced
- `inlabel::String=""`: label of files used for initial conditions
- `fittedparam::Vector=Int[]`: vector of rate indices to be fit, e.g. [1,2,3,5,6,7]
- `fixedeffects::Tuple=tuple()`: tuple of vectors of rates that are fixed where first index is fit and others are fixed to first, e.g. ([3,8],) means index 8 is fixed to index 3
     (only first parameter should be included in fittedparam)
- `transitions::Tuple=([1,2],[2,1])`: tuple of vectors that specify state transitions for G states, e.g. ([1,2],[2,1]) for classic 2 state telegraph model and ([1,2],[2,1],[2,3],[3,1]) for 3 state kinetic proof reading model
- `G::Int=2`: number of gene states
- `R::Int=0`: number of pre-RNA steps (set to 0 for classic telegraph models)
- `S::Int=0`: number of splice sites (set to 0 for classic telegraph models and R for GRS models, when insertstep > 1, S is actually R - insertstep + 1)
- `insertstep::Int=1`: R step where reporter is inserted
- `Gfamily=""`: String describing type of G transition model, e.g. "3state", "KP" (kinetic proofreading), "Refractory"
- `root="."`: name of root directory for project, e.g. "scRNA"
- `priormean=Float64[]`: mean of prior rate distribution
- 'priorcv=10.`: coefficient of variation for the rate prior distributions, default is 10.
- `nalleles=2`: number of alleles, value in alleles folder will be used if it exists  
- `onstates=Int[]`: vector of on or sojourn states
- `decayrate=1.0`: decay rate of mRNA, if set to -1, value in halflives folder will be used if it exists
- `splicetype=""`: RNA pathway for GRS models, (e.g. "offeject" =  spliced intron is not viable)
- `probfn=prob_GaussianMixture`: probability function for hmm observation probability (e.g. prob_GaussianMixture)
- `noiseparams=5`: number of parameters of probfn
- `weightind=5`: parameter index of bias probability of mixtures, e.g. noiseparams=5, weightind=5 means last noise parameter is for mixture bias
- `hierarchical=tuple()`: tuple of hierchical model parameters
- `ratetype="median"`: which rate to use for initial condition, choices are "ml", "mean", "median", or "last"
- `propcv=0.01`: coefficient of variation (mean/std) of proposal distribution, if cv <= 0. then cv from previous run will be used
- `maxtime=Float64=60.`: maximum wall time for run, default = 60 min
- `samplesteps::Int=1000000`: number of MCMC sampling steps
- `warmupsteps=0`: number of MCMC warmup steps to find proposal distribution covariance
- `annealsteps=0`: number of annealing steps (during annealing temperature is dropped from tempanneal to temp)
- `temp=1.0`: MCMC temperature
- `tempanneal=100.`: annealing temperature
- `temprna=1.`: reduce rna counts by temprna compared to dwell times
- `burst=false`: if true then compute burst frequency
- `optimize=false`: use optimizer to compute maximum likelihood value
- `writesamples=false`: write out MH samples if true, default is false

```

```
    makeswarm(genes::Vector{String}; <keyword arguments> )


write swarm and fit files used on biowulf
creates a run for each gene

#Arguments
- `genes`: vector of genes
- `nchains::Int=2`: number of MCMC chains = number of processors called by Julia, default = 2
- 'nthreads::Int=1`: number of Julia threads per processesor, default = 1
- `swarmfile::String="fit"`: name of swarmfile to be executed by swarm
- `batchsize=1000`: number of jobs per swarmfile, default = 1000
- `juliafile::String="fitscript`: name of file to be called by julia in swarmfile
- `thresholdlow::Float=0`: lower threshold for halflife for genes to be fit
- `threhsoldhigh::=Inf`: upper threshold
- `datatype::String=""`: String that desecribes data type, choices are "rna", "rnaonoff", "rnadwelltime", "trace", "tracenascent", "tracerna"
- `dttype=String[]`: types are "OFF", "ON", for R states and "OFFG", "ONG" for G states
- `datapath=""`: path to data file or folder or array of files or folders
- `cell::String=""': cell type for halflives and allele numbers
- `datacond=""`: string or vector of strings describing data treatment condition, e.g. "WT", "DMSO" or ["DMSO","AUXIN"]
- `traceinfo=(1., 1., .65)`: vector of frame traceinfo of intensity traces and transient time
- `nascent=(1,2,1.)`: vector of number of spots, and total number of locations (e.g. number of cells times number of alleles/cell)
- `infolder::String=""`: result folder used for initial parameters
- `resultfolder::String=test`: folder for results of MCMC run
- `label::String=""`: label of output files produced
- `inlabel::String=""`: label of files used for initial conditions
- `fittedparam::Vector=Int[]`: vector of rate indices to be fit, e.g. [1,2,3,5,6,7]
- `fixedeffects=tuple()`: tuple of vectors of rates that are fixed where first index is fit and others are fixed to first, e.g. ([3,8],) means  index 8 is fixed to index 3
     (only first parameter should be included in fixedeffects)
- `transitions::Tuple=([1,2],[2,1])`: tuple of vectors that specify state transitions for G states, e.g. ([1,2],[2,1]) for classic 2 state telegraph model and ([1,2],[2,1],[2,3],[3,1]) for 3 state kinetic proof reading model
- `G::Int=2`: number of gene states
- `R::Int=0`: number of pre-RNA steps (set to 0 for classic telegraph models)
- `S::Int=0`: number of splice sites (set to 0 for classic telegraph models and R for GRS models, when insertstep > 1, S is R - insertstep + 1))
- `insertstep::Int=1`: R step where reporter is inserted
- `Gfamily=""`: String describing type of G transition model, e.g. "3state", "KP" (kinetic proofreading), "Refractory"
- `root="."`: name of root directory for project, e.g. "scRNA"
- `priormean=Float64[]`: mean of prior rate distribution
- 'priorcv=10.`: coefficient of variation for the rate prior distributions, default is 10.
- `nalleles=2`: number of alleles, value in alleles folder will be used if it exists  
- `onstates=Int[]`: vector of on or sojourn states, e.g. [[2,3],[]], use empty vector for R states, do not use [] for R=0 models
- `decayrate=1.0`: decay rate of mRNA, if set to -1, value in halflives folder will be used if it exists
- `splicetype=""`: RNA pathway for GRS models, (e.g. "offeject" =  spliced intron is not viable)
- `probfn=prob_GaussianMixture`: probability function for hmm observation probability (e.g. prob_GaussianMixture)
- `noiseparams=5`: number of parameters of probfn
- `weightind=5`: parameter index of bias probability of mixtures, e.g. noiseparams=5, weightind=5 means last noise parameter is for mixture bias
- `hierarchical=tuple()`: tuple of hierchical model parameters
- `ratetype="median"`: which rate to use for initial condition, choices are "ml", "mean", "median", or "last"
- `propcv=0.01`: coefficient of variation (mean/std) of proposal distribution, if cv <= 0. then cv from previous run will be used
- `maxtime=Float64=60.`: maximum wall time for run, default = 60 min
- `samplesteps::Int=1000000`: number of MCMC sampling steps
- `warmupsteps=0`: number of MCMC warmup steps to find proposal distribution covariance
- `annealsteps=0`: number of annealing steps (during annealing temperature is dropped from tempanneal to temp)
- `temp=1.0`: MCMC temperature
- `tempanneal=100.`: annealing temperature
- `temprna=1.`: reduce rna counts by temprna compared to dwell times
- `burst=false`: if true then compute burst frequency
- `optimize=false`: use optimizer to compute maximum likelihood value
- `writesamples=false`: write out MH samples if true, default is false
- `src=""`: path to folder containing StochasticGene.jl/src
```




```
run_mh(data,model,options)

returns fits, stats, measures

Run Metropolis-Hastings MCMC algorithm and compute statistics of results

-`data`: AbstractExperimentalData structure
-`model`: AbstractGmodel structure with a logprior function
-`options`: MHOptions structure

model and data must have a likelihoodfn function
```
```
    write_dataframes(resultfolder::String, datapath::String; measure::Symbol=:AIC, assemble::Bool=true, fittedparams=Int[])

  write_dataframes(resultfolder::String,datapath::String;measure::Symbol=:AIC,assemble::Bool=true)

  collates run results into a csv file

Arguments
- `resultfolder`: name of folder with result files
- `datapath`: name of folder where data is stored
- `measure`: measure used to assess winner
- `assemble`: if true then assemble results into summary files
```
```
    simulator(r::Vector{Float64}, transitions::Tuple, G::Int, R::Int, S::Int, insertstep::Int; nalleles::Int=1, nhist::Int=20, onstates::Vector=Int[], bins::Vector=Float64[], traceinterval::Float64=0.0, probfn=prob_GaussianMixture, noiseparams::Int=5, totalsteps::Int=1000000000, totaltime::Float64=0.0, tol::Float64=1e-6, reporterfn=sum, splicetype="", verbose::Bool=false)

Simulate any GRSM model. Returns steady state mRNA histogram. If bins not a null vector will return a vector of the mRNA histogram and ON and OFF time histograms. If traceinterval > 0, it will return a vector containing the mRNA histogram and the traces

#Arguments
- `r`: vector of rates
- `transitions`: tuple of vectors that specify state transitions for G states, e.g. ([1,2],[2,1]) for classic 2 state telegraph model and ([1,2],[2,1],[2,3],[3,1]) for 3 state kinetic proof reading model
- `G`: number of gene states
- `R`: number of pre-RNA steps (set to 0 for classic telegraph models)
- `S`: number of splice sites (set to 0 for G (classic telegraph) and GR models and R for GRS models)
- `insertstep`: reporter insertion step
 	
#Named arguments
- `nalleles`: Number of alleles
- `nhist::Int`: Size of mRNA histogram
- `onstates::Vector`: a vector of ON states (use empty set for any R step is ON) or vector of vector of ON states
- `bins::Vector=Float64[]`: vector of time bins for ON and OFF histograms or vector of vectors of time bins
- `probfn`=prob_GaussianMixture: reporter distribution
- `traceinterval`: Interval in minutes between frames for intensity traces.  If 0, traces are not made.
- `totalsteps`::Int=10000000: maximum number of simulation steps (not usred when simulating traces)
- `tol`::Float64=1e-6: convergence error tolerance for mRNA histogram (not used when simulating traces are made)
- `totaltime`::Float64=0.0: total time of simulation
- `splicetype`::String: splice action
- `reporterfn`=sum: how individual reporters are combined
- `verbose::Bool=false`: flag for printing state information
    
#Example:

julia> h=simulator(r,transitions,3,2,2,1,nhist=150,bins=[collect(5/3:5/3:200),collect(.1:.1:20)],onstates=[Int[],[2,3]],nalleles=2)

```
```
    simulate_trace_data(datafolder::String;ntrials::Int=10,r=[0.038, 1.0, 0.23, 0.02, 0.25, 0.17, 0.02, 0.06, 0.02, 0.000231,30,20,200,100,.8], transitions=([1, 2], [2, 1], [2, 3], [3, 1]), G=3, R=2, S=2, insertstep=1,onstates=Int[], interval=1.0, totaltime=1000.)

create simulated trace files in datafolder
```

```
    predicted_trace(data::Union{AbstractTraceData,AbstractTraceHistogramData}, model)

return predicted traces of fits using Viterbi algorithm
```