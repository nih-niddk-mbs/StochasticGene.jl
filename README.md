# StochasticGene.jl

Julia package to simulate and fit stochastic models of gene transcription to experimental data. The data acceptable include distributions of mRNA counts per cell (e.g. single molecule FISH (smFISH) or single cell RNA sequence (scRNA) data)), image intensity traces from live cell imaging, dwell time distributions of reporters (e.g. MS2 or PP7) imaged in live cells, and in combination. The transcription models considered are stochastic (continuous Markov systems) with an arbitrary number of G (gene) states, R (pre-RNA) steps, S splice sites (up to R), and reporter insertion step (insertstep). The gene always occupies one of the G states and there are random (user specified) transitions between G states.  One of the G states is an active state where transcription can be initiated and the first R step becomes occupied. An irreversible forward transition can then occur to the next R step if that step is unoccupied mimicking elongation. An mRNA molecule is ejected from the final (termination) R step where it then decays stochastically. The model can account for multiple alleles of the gene in the same cell. The user can specify which R step is considered visible when occupied (i.e. contains a reporter). For example, if the pre-RNA MS2 construct is inserted into an intron that is transcribed early then you could make the reporter insertstep to be 1 and all R steps will be visible. But if the reporter is transcribed late then you could choose a later R step, like step 3, and only R steps at and after step 3 are considered visible. The genes/alleles can also be coupled. In the original model in Rodriguez et al. Cell (2018), the reporter was in an exon and was visible at all steps and then ejected. In Wan et al. Cell (2021), the reporter is inserted into an intron and thus can be spliced out before the polymerase reaches the final R step. Models are allowed to have no R steps (i.e. classic telegraph models but with arbitrary numbers of G states (rather thabn usual two)) where an mRNA molecule can be stochastically produced when the gene occupies the active G state. You can also specify reporters for G states (e.g. transcription factors) and more than one simultaneous reporter (e.g. reporters for introns and transcription factors). The model can also be run on multiple processors and has been used to fit scRNA from the entire genome. Trzaskoma et al. Science Advances (2024).

The package has functions to specify the models, prepare the data, compute the predicted data, apply a Metropolis-Hastings markov chain monte carlo (MCMC) algorithm to fit the parameters of the models to the data (i.e. compute Bayesian posterior distributions), explicitly simulate models, and analyze the results. Unfortunately, not all capabilities are documented so just send me a message if you have any questions.

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

After the package has installed, you can check if it works by running (you can skip this part if you don't want to wait)

```
julia> ] test StochasticGene
```
(this could take 20 or more minutes)
Currently, test may throw some warnings or errors regarding dependencies. These are from external packages that StochasticGene references but does not rely on so it can be ignored.

StochasticGene will be updated on Github periodically. To update to a new version type

```
julia> ] update StochasticGene
```

To install StochasticGene on Biowulf, log on and at the prompt type:

```
[username@biowulf ~]$ sinteractive --constraint=x2695 --mem=64G
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

Also, Julia is periodically updated on Biowulf and StochasticGene will not carry forward to the latest version.  You can still call previous versions when you allocated CPUs or reinstall StochasticGene in the new version.

### Creating folder structure to run StochasticGene

StochasticGene requires a specific directory or folder structure where data are stored and results are saved.  At the top is the `root` folder (e.g. "scRNA" or "RNAfits") with subfolders `data` and `results`. Inside `data` are two more folders  `alleles` and `halflives`,  containing allele numbers and half lives, respectively.  The command 

```
julia> using StochasticGene

julia> rna_setup("scRNA")
```
will create the folder structure in the folder `scRNA` with the `data` and `results` subfolders, with some mock data in the `data` folder. Typing `rna_setup()` will assume the current folder is the root. Files for allele numbers and halflives for desired cell types can be added directly to the `data` folder.  These files should be csv format and have the form `[cell name]_alleles.csv` or `[cell name]_halflife.csv`.  Halflives for the then genes in Wan et al. (2021) for HBEC cells are built into the system.

To move into the folder type

```
julia> cd("scRNA")
```
To confirm you are in the correct directory, use
```
julia> pwd()
```

### Fitting data with StochasticGene

To fit a model, you need to load data and choose a model. Data currently accepted are stationary histograms (e.g. smFISH or scRNA), intensity time traces (e.g. trk files), and dwell time distributions. Different data types from the same experiment can also be fit simultaneously. For example, RNA histograms from smFISH can be fit together with multiple intensity traces or dwell time histograms from live cell recordings.

Models are distringuished by the number of G states, the transition graph between G states (i.e. which G states can transitions to which other G states), number of R steps, the step where the reporter is inserted, and whether splicing occurs.  For intensity traces and dwell time distributions, the sojourn or "on" states (i.e. R steps or G states where the reporter is visible) must also be specified. Multiple sets of on states are allowed.

Data are fit using the `fit` function, which has named arguments (see API below) to determine the type of data to be fit, where it can be found, what model to be used, and what options to be set for the fitting algorithm. The default, run by typing `fit()`, will fit the mock rna histogram data installed by `rna_setup(root)` with a simple two state telegraph model (2 G stsates, no R steps).  Data is in the folder "root/data/HCT116_testdata", where root is specified by the user (not the root of the computer system but the folder where the data and results subfolders reside). The default root is the folder in which julia was launched but can be specified using the named argument `root`. The `fit` function returns six variables. This will only work if you start Julia in the root folder.

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

During the run, the `fit` function prints out some of the information in fits, stats, and measures structures. `fit` will look for previous results in the folder specified by the argument `infolder` as an initial condition to start the MCMC run. In this case, there were no previous runs and thus the default was used, which is the priors of the parameters. In this run three rates were fit (set by argument `fittedparams`); the median of the posterior and the maximum likelihood parameters of the resulting fits are printed. The run was capped to a maximum of 60 seconds of real time (set by the argument `maxtime`) and `Acceptance` shows that out of 613219 MCMC samples, 325538 were accepted by the MCMC algorithm. The positive number `Deviance` is the difference in likelihood between a perfect fit and the resulting fit (for histograms only). `rhat` is a measure of MCMC convergence. It is always larger than 1 with 1 being ideal. 

The two state telegraph model has 4 transition rates, which are stored in a single vector. The order is 1) G state transitions (in the order as specified by the transition vector), 2) eject rate, 3) decay rate.  For general GRSM models, the order is 1) G state transitions 2) R transitions, 3) S transitions, 4) decay.  If there is no splicing then the S transitions are left out.  Not all the rates need to be fitted. In the above example, the decay rate is not fitted. This is specified by the fittedparams argument (see API). The posterior median, mean, standard deviations, and MADs will only be for the fitted params.

In this particular run, only one chain was used so the measures are not very informative. To use more chains, specify more processors with

```
bash> julia -p 4

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

The `datatype` argument is a String that specifies the types of data to be fit. Currently, there are six choices: 1) "rna", which expects a single rna histogram file where the first column is the histogram; 2) "rnaonoff", which expects a rna histogram file and a three column dwelltime file with columns: bins, ON time distribution, and OFF time distribution; 3) "rnadwelltime", which fits an rna histogram together with multiple dwell time histograms specified by a vector of dwell time types, the choices being "ON","OFF" for R step reporters and "ONG", "OFFG" for G state reporters. Each dwelltime histogram is a two column file of the bins and dwell times and datapath is a vector of the paths to each; 4) "trace", which expects a folder of trace intensity (e.g. trk) files; 5) "tracerna": trace files with RNA histogram; 6) "tracejoint", which fits two simultaneous recorded traces using a coupled model. The data files are specified by the argument `datapath`. This can be a string pointing to a file or a folder containing the file.  If it points to a folder then `fit` will use the arguments `gene` and `datacond` to identify the correct file. For datatypes that require more than one data file, `datapath` is a vector of paths.  The path can include or not include the root folder.  

### Example fitting traces

You will still need the same folder structure with a root folder with results and data subfolders.

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

which generates 10 mock trk files in folder "data/testtraces". See API below for the parameters used in the simulated data.  We fit a model with a G=3, R=2, S=2, insertstep=1 model using the command below. This model assumes that there are 3 Gene states and 2 pre-RNA steps with the reporter inserted at the first step. The reporter can also be spliced out at any step. If S < R then splicing will only occur in the first S available slots, starting from 1. Insertstep attempts to account for the fact that the reporter could be inserted anywhere along the gene (e.g. 3' or 5' end). If the signal intensity indicates more than one reporter is present during a burst then choose insertstep to be smaller than R. However, if bursts seem to only contain one reporter then use insertstep = R. The `transitions` keyword takes a tuple of vectors where each vector indicates a G transition to be included. For example [1, 2] means that there exists a transition from G state 1 to G state 2.

```

julia> fits, stats, measures, data, model, options = fit(datatype="trace",nchains=4,transitions=([1, 2], [2, 1], [2, 3], [3, 1]), G=3, R=2, S=2, insertstep=1, datapath="data/testtraces", gene="test", datacond="", decayrate=1.,noisepriors=[40., 20., 200., 10.]);
HCT116_test not found
./results/HCT116_test/rates_trace-HCT116-nstate__test_3221_1.txt does not exist
2025-01-10T11:24:23.788
Gene: testLabel:  G R S insertstep: 3221
data: data/testtraces
in: HCT116_test out: HCT116_test
maxtime: 60.0
10
data/testtraces

(1.0, 1, -1, 1.0)
No rate file, set rate to prior
[0.01, 0.01, 0.01, 0.01, 0.1, 0.3333333333333333, 0.3333333333333333, 0.1, 0.1, 1.0, 40.0, 20.0, 200.0, 10.0]
initial ll: 80283.84924591208
final max ll: 52072.877782295
median ll: 59649.5738579672
Median fitted rates: [0.010990214988614367, 0.007759982708079304, 0.009675101244773334, 0.00647827439174417, 0.35238021798551666, 1.8043372460110718, 1.7175355401214196, 0.1628658491792051, 0.10206281750564142, 32.161628025595455, 51.33142446945536, 209.43545981757265, 14.817798527873652]
ML rates: [0.010501832166699481, 0.005916260295568495, 0.011746187235209794, 0.00465431537381948, 1.0946492144589024, 4.509832072463753, 2.2247620219458835, 0.20626061166543908, 0.10502518236340562, 3.011698651381525, 4.642095991085645, 156.64978288985006, 96.59297507554737]
Acceptance: 4424/9279
rhat: 2.541388750409241
2025-01-10T11:25:42.925


```
The results will be found in a folder called results/HCT116_test, which you can set using the `resultfolder` keyword (see the API for fit below). Each time you run fit, it will first check 'infolder' for a "rate" file to load in the results of previous fits. Thus, setting `infolder` and `resultfolder` to the same folder will allow you to iterate the fits. When fitting traces, it is imperative that you give very good priors for the background and signal mean intensities.  These are specified by the `noisepriors` keyword. The default is to fit both the background and the signal intensities with a Gaussian distribution and thus the 4 parameters in `noisepriors` corresponds to the [mean of the background, S.D. of the background, mean of the signal, S.D. of the signal]. The default distribution assumes that the signal intensities are additive so if two reporters are present, the the signal distribution will have a mean twice as large and variance twice as large. Other distributions can be selected using the `probfn` keyword. Options include nonadditive Gaussians and mixed Gaussians. Some will require more parameters. Generally, the default will be sufficient. It is also sometimes beneficial to not fit the signal intensities. Which parameters are to be fit can be specified by the `fittedparam` keyword.

In this example run, `rhat` is above 2.5 indicating that the number of samples was probably insufficient to obtain a good sampling of the posterior distributions. Either `maxtime` or `samplesteps` needs to be increased.


### Batch fitting on Biowulf using `swarm`.

(You made need to run a few times on Biowulf after installing or updating StochasticGene before it works as it sometimes takes too long to precompile on the first use)

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
`-t` flag is for time
`-g` flag is for memory

(choose a time longer than maxtime (remember to convert seconds to hours), and make sure to allocate enough memory)

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

For coupled transcribing units, arguments transitions, G, R, S, insertstep, and trace become tuples of the single unit type, e.g. If two types of transcription models are desired with G= 2 and G=3 then
then G = (2,3). 
#Arguments
- `annealsteps=0`: number of annealing steps (during annealing temperature is dropped from tempanneal to temp)
- `burst=false`: if true then compute burst frequency
- `cell::String=""`: cell type for halflives and allele numbers
- `coupling=tuple()`: if nonempty, a 4-tuple where elements are 
    1. tuple of model indices corresponding to each unit, e.g. (1, 1, 2) means that unit 1 and 2 use model 1 and unit 3 uses model 2
    2. tuple of vectors indicating source units for each unit, e.g. ([2,3], [1], Int[]) means unit 1 is influenced by source units 2 and 3, unit 2 is influenced by unit 1 and unit 3 is uninfluenced.
    3. source states: tuple of vectors of strings, e.g. (["G3","R1"], []) means that model 1 influences other units whenever it is in G state 3 or R step 1 (if a number is not included (e.g. (R,0)) then all the Gstates or R steps are included), 
        while model 2 does not influence any other unit
    4. target transitions:tuple, e.g. ([], 4) means that model 1 is not influenced by any source while model 2 is influenced by sources at transition 4. Transitions
        are number consecutively by order of the transition rates. So for a G=2 model, transition 1 is the G1 to G2 transition and transition 3 is the initiation transition
    5. Int indicating number of coupling parameters
- `datatype::String=""`: String that describes data type, choices are "rna", "rnaonoff", "rnadwelltime", "trace", "tracerna", "tracejoint", "tracegrid"
- `datacond=""`: string or vector of strings describing data, e.g. "WT", "DMSO" or ["DMSO","AUXIN"], ["gene","enhancer"]
- `datapath=""`: path to data file or folder or array of files or folders
- `decayrate=1.0`: decay rate of mRNA, if set to -1, value in halflives folder will be used if it exists
- `dttype=String[]`: dwelltime types, choices are "OFF", "ON", for R states and "OFFG", "ONG" for G states
- `elongationtime=6.0`: average time for elongation, vector of times for coupled model
- `fittedparam::Vector=Int[]`: vector of rate indices to be fit, e.g. [1,2,3,5,6,7]  (applies to shared rates for hierarchical models, fitted hyper parameters are specified by individual fittedparams)
- `fixedeffects::String`: if "fixed" is included after a hyphen, then fixedeffects Tuple will be created such that R transitions are fixed to be identical
- `fixedeffects::Tuple=tuple()`: tuple of vectors of rates that are fixed where first index is fit and others are fixed to first, e.g. ([3,8],) means index 8 is fixed to index 3 (only first parameter should be included in fittedparam) (applies to shared rates for hierarchical models)
- `gene::String="MYC"`: gene name
- `grid=nothing`: Int number of grid points for grid model
- `G=2`: number of gene states, for coupled models G, R, S, and insertstep are vectors (vector for coupled models)
- `hierarchical=tuple()`: empty tuple for nonhierarchical model; 3-tuple for hierarchical: hierarchical=(number of hyper parameter sets::Int, individual fittedparams::Vector, individual fixedeffects::Tuple),
    for hierarchical models the keywords `fittedparam` and `fixedeffects` pertain to shared rates.  rates are given by a single vector that can be reshaped into a matrix where the columns correspond to the model rates and noise params, the first nhyper rows pertain to the shared and hyper parameter rates (whether fit or not), 
    usually the first row is the shared and mean hyper parameters and the 2nd are the standard deviations, the rest of the rows are the individual rates and noise params
- `infolder::String=""`: result folder used for initial parameters
- `inlabel::String=""`: label of files used for initial conditions
- `insertstep=1`: R step where reporter is inserted
- `label::String=""`: label of output files produced
- `maxtime=Float64=60.`: maximum wall time for run, default = 60 min
- `method=Tsit5()`: DifferentialEquations.jl numerical method (e.g. Tsit5(), lsoda(),...); use a tuple for hierarchical models: method = tuple(method, Bool) = (numerical method (currently not used), true if transition rates are shared)
- `nalleles=1`: number of alleles, value in alleles folder will be used if it exists, always set to 1 if coupling nonempty (multiple alleles handled as extra identical units)
- `nchains::Int=2`: number of MCMC chains = number of processors called by Julia, default = 2
- `noisepriors=[]`: priors of observation noise (use empty set if not fitting traces), superseded if priormean is set
- `onstates::Vector{Int}=Int[]`: vector of on or sojourn states, e.g. [[2,3],Int[]], use empty vector for R states, do not use Int[] for R=0 models
- `optimize=false`: use optimizer to compute maximum likelihood value
- `priormean=Float64[]`: mean rates of prior distribution (must set priors for all rates including those that are not fitted)
- `priorcv=10.`: (vector or number) coefficient of variation(s) for the rate prior distributions, default is 10.
- `probfn=prob_Gaussian`: probability function for HMM observation probability (i.e., noise distribution)
- `propcv=0.01`: coefficient of variation (mean/std) of proposal distribution, if cv <= 0. then cv from previous run will be used
- `resultfolder::String=test`: folder for results of MCMC run
- `R=0`: number of pre-RNA steps (set to 0 for classic telegraph models)
- `root="."`: name of root directory for project, e.g. "scRNA"
- `samplesteps::Int=1000000`: number of MCMC sampling steps
- `S=0`: number of splice sites (set to 0 for classic telegraph models and R - insertstep + 1 for GRS models)
- `splicetype=""`: RNA pathway for GRS models, (e.g., "offeject" = spliced intron is not viable)
- `temp=1.0`: MCMC temperature
- `tempanneal=100.`: annealing temperature
- `temprna=1.`: reduce RNA counts by temprna compared to dwell times
- `traceinfo=(1.0, 1., -1, 1.)`: 4-tuple = (frame interval of intensity traces, starting frame time in minutes, ending frame time (use -1 for last index), fraction of observed active traces); for simultaneous joint traces, the fraction of active traces is a vector of the active fractions for each trace, e.g. (1.0, 1., -1, [.5, .7]) 
- `TransitionType=""`: String describing G transition type, e.g. "3state", "KP" (kinetic proofreading), "cyclic", or if hierarchical, coupled
- `transitions::Tuple=([1,2],[2,1])`: tuple of vectors that specify state transitions for G states, e.g. ([1,2],[2,1]) for classic 2-state telegraph model and ([1,2],[2,1],[2,3],[3,1]) for 3-state kinetic proofreading model
- `warmupsteps=0`: number of MCMC warmup steps to find proposal distribution covariance
- `writesamples=false`: write out MH samples if true, default is false

Example:

If you are in the folder where data/HCT116_testdata is installed, then you can fit the mock RNA histogram running 4 mcmc chains with

bash> julia -p 4

julia> fits, stats, measures, data, model, options = fit(nchains = 4)

```


```
    makeswarm(;<keyword arguments>)

write swarm and fit files used on biowulf


#Arguments

- 'nthreads=1`: number of Julia threads per processesor, default = 1
- `swarmfile::String="fit"`: name of swarmfile to be executed by swarm
- `juliafile::String="fitscript`: name of file to be called by julia in swarmfile
- `src=""`: path to folder containing StochasticGene.jl/src (only necessary if StochasticGene not installed)

and all keyword arguments of function fit(; <keyword arguments> )

see fit

Note:: the keyword 'method' is handled slightly differently here than in the function fit.  In fit it is a function (i.e. numerical method function used in DifferentialEquations.jl) or a tuple
of the numerical method and a Bool for hierarchical models.  However, in biowulf.jl, the numerical method must be a String, i.e. use "lsoda()" for lsoda().  This is because when Julia writes
the function, it will parse the numerical method rather than just writing it.

```

```
    makeswarm_genes(genes::Vector{String}; <keyword arguments> )

write a swarmfile and fit files to run all each gene in vector genes

# Arguments
- `genes`: vector of genes
- `batchsize=1000`: number of jobs per swarmfile, default = 1000

and all arguments in makeswarm(;<keyword arguments>)


    Examples

julia> genes = ["MYC","SOX9"]

julia> makeswarm(genes,cell="HBEC")

```

```
    makeswarm_genes(;<keyword arguments> )

@JuliaRegistrator register()
 
#Arguments
    - `thresholdlow::Float=0`: lower threshold for halflife for genes to be fit
    - `threhsoldhigh::=Inf`: upper threshold

    and all keyword arguments in makeswarm(;<keyword arguments>)
```

```
    makeswarm(models::Vector{ModelArgs}; <keyword arguments> )

creates a run for each model

#Arguments
- `models::Vector{ModelArgs}`: Vector of ModelArgs structures

and all keyword arguments in makeswarm(;<keyword arguments>)

```

```
run_mh(data,model,options)

returns fits, stats, measures

Run Metropolis-Hastings MCMC algorithm and compute statistics of results

-`data`: AbstractExperimentalData structure
-`model`: AbstractGeneTransitionModel structure with a logprior function
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
   simulator(r, transitions, G, R, S, insertstep; coupling=tuple(), nalleles=1, nhist=20, onstates=Int[], bins=Float64[], traceinterval::Float64=0.0, probfn=prob_Gaussian, noiseparams=4, totalsteps::Int=10000, totaltime::Float64=0.0, tol::Float64=1e-6, reporterfn=sum, splicetype="", verbose::Bool=false)

Simulate any GRSM model. Returns steady state mRNA histogram. If bins not a null vector will return a vector of the mRNA histogram and ON and OFF time histograms. If traceinterval > 0, it will return a vector containing the mRNA histogram and the traces

#Arguments
- `r`: vector of rates
- `transitions`: tuple of vectors that specify state transitions for G states, e.g. ([1,2],[2,1]) for classic 2 state telegraph model and ([1,2],[2,1],[2,3],[3,1]) for 3 state kinetic proof reading model
- `G`: number of gene states
- `R`: number of pre-RNA steps (set to 0 for classic telegraph models)
- `S`: number of splice sites (set to 0 for G (classic telegraph) and GR models and R for GRS models)
- `insertstep`: reporter insertion step
 	
#Named arguments
- `bins::Vector=Float64[]`: vector of time bin vectors for each set of ON and OFF histograms or vector of vectors of time bins (one time bin vector for each onstate)
- `coupling=tuple()`: if nonempty, a 4-tuple where elements are 
    1. tuple of model indices corresponding to each unit, e.g. (1, 1, 2) means that unit 1 and 2 use model 1 and unit 3 uses model 2
    2. tuple of vectors indicating source units for each unit, e.g. ([2,3], [1], Int[]) means unit 1 is influenced by source units 2 and 3, unit 2 is influenced by unit 1 and unit 3 is uninfluenced.
    3. source states, e.g. (3,0) means that model 1 influences other units whenever it is in G state 3, while model 2 does not influence any other unit
    4. target transitions, e.g. (0, 4) means that model 1 is not influenced by any source while model 2 is influenced by sources at G transition 4.
    5. Int indicating number of coupling parameters
- `nalleles`: Number of alleles, set to 1 if coupling nonempty
- `nhist::Int`: Size of mRNA histogram
- `onstates::Vector`: a vector of vector of ON states (use empty set for any R step is ON), ON and OFF time distributions are computed for each ON state set
- `probfn`=prob_Gaussian: reporter distribution
- `reporterfn`=sum: how individual reporters are combined
- `splicetype`::String: splice action
- `tol`::Float64=1e-6: convergence error tolerance for mRNA histogram (not used when simulating traces are made)
- `totalsteps`::Int=10000000: maximum number of simulation steps (not usred when simulating traces)
- `totaltime`::Float64=0.0: total time of simulation
- `traceinterval`: Interval in minutes between frames for intensity traces.  If 0, traces are not made.
- `verbose::Bool=false`: flag for printing state information
    
#Example:

julia> h=simulator([.1, .1, .1, .1, .1, .1, .1, .1, .1, .01],([1,2],[2,1],[2,3],[3,2]),3,2,2,1,nhist=20,bins=[collect(2.:2.:200),collect(.2:.2:20)],onstates=[Int[],[3]],nalleles=2, totalsteps = 200000)
5-element Vector{Vector}:
 [7823.967526508377, 33289.19787176562, 69902.6774554014, 92942.59127561412, 91816.91189325438, 70259.88319069796, 43895.28637479579, 22426.725922619895, 9005.190755732247, 3043.2417332890695, 1005.6412773072143, 203.11430396725336, 6.420639815427421, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
 [2911, 2568, 2228, 1694, 1354, 1088, 819, 661, 514, 401  …  0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
 [622, 643, 592, 596, 553, 524, 520, 489, 448, 437  …  19, 13, 12, 12, 14, 10, 10, 17, 8, 8]
 [550, 584, 569, 555, 498, 510, 497, 495, 489, 487  …  89, 96, 107, 99, 89, 103, 86, 97, 87, 77]
 [593, 519, 560, 512, 492, 475, 453, 468, 383, 429  …  84, 73, 85, 92, 73, 81, 85, 101, 79, 78]
 
```
```
    simulate_trace_data(datafolder::String;ntrials::Int=10,r=[0.038, 1.0, 0.23, 0.02, 0.25, 0.17, 0.02, 0.06, 0.02, 0.000231,30,20,200,100,.8], transitions=([1, 2], [2, 1], [2, 3], [3, 1]), G=3, R=2, S=2, insertstep=1,onstates=Int[], interval=1.0, totaltime=1000.)

create simulated trace files in datafolder
```
