# This file is part of StochasticGene.jl
#
#genetrap.jl
#Fit GRSM models to live cell and smFISH data


"""
Gene information
"""
genes_gt() = ["CANX"; "DNAJC5"; "ERRFI1"; "KPNB1"; "MYH9"; "Rab7a"; "RHOA"; "RPAP3"; "Sec16A"; "SLC2A1"]
genelength_gt() = Dict([("Sec16A", 42960); ("SLC2A1", 33802); ("ERRFI1", 14615); ("RHOA", 52948); ("KPNB1", 33730); ("MYH9", 106741); ("DNAJC5", 40930); ("CANX", 32710); ("Rab7a", 88663); ("RPAP3", 44130)])
MS2end_gt() = Dict([("Sec16A", 5220); ("SLC2A1", 26001); ("ERRFI1", 5324); ("RHOA", 51109); ("KPNB1", 24000); ("MYH9", 71998); ("DNAJC5", 14857); ("CANX", 4861); ("Rab7a", 83257); ("RPAP3", 38610)])
halflife_gt() = Dict([("CANX", 50.0), ("DNAJC5", 5.0), ("ERRFI1", 1.35), ("KPNB1", 9.0), ("MYH9", 10.0), ("Rab7a", 50.0), ("RHOA", 50.0), ("RPAP3", 7.5), ("Sec16A", 8.0), ("SLC2A1", 5.0)])


"""
    fit_genetrap()

"""

function fit_genetrap(nchains, maxtime, gene::String, transitions, G::Int, R::Int, S::Int, insertstep::Int; onstates=[], priorcv=10.0, propcv=0.01, fittedparam=collect(1:num_rates(transitions, R, R, insertstep)-1), infolder::String="test", folder::String="test", samplesteps::Int=1000, nalleles::Int=2, label="gt", rnatype="", warmupsteps=0, annealsteps=0, temp=1.0, tempanneal=100.0, tempfish=1.0, root::String=".", burst=false)
    println(now())
    folder = folder_path(folder, root, "results", make=true)
    infolder = folder_path(infolder, root, "results")
    data, model = genetrap(root, gene, transitions, G, R, S, insertstep, 2, rnatype, fittedparam, infolder, folder, label, "ml", tempfish, priorcv, propcv, onstates)
    println("size of histogram: ", data.nRNA)
    options = MHOptions(samplesteps, warmupsteps, annealsteps, maxtime, temp, tempanneal)
    println(model.rates)
    print_ll(data, model)
    fits, stats, measures = run_mh(data, model, options, nchains)
    if burst
        bs = burstsize(fit, model)
    else
        bs = 0
    end
    # finalize(data,model,fits,stats,measures,1.,folder,0,bs,root)
    finalize(data, model, fits, stats, measures, temp, folder, 0, burst, false, root)
    println(now())
    return data, model_genetrap(gene, get_rates(fits.parml, model), transitions, G, R, S, insertstep, fittedparam, nalleles, data.nRNA + 2, priorcv, propcv, onstates, rnatype), fits, stats, measures
end

"""
genetrap()

Load data, model and option structures for metropolis-hastings fit
root is the folder containing data and results

FISH counts is divided by tempfish to adjust relative weights of data
set tempfish = 0 to equalize FISH and live cell counts, tempfish < 0 for FISH alone
"""
function genetrap(root, gene::String, transitions::Tuple, G::Int, R::Int, S::Int, insertstep::Int, nalleles, rnatype::String, fittedparam::Vector, infolder::String, resultfolder::String, label::String, rtype::String, tempfish, priorcv, propcv, onstates)
    r = readrates_genetrap(infolder, rtype, gene, label, G, R, S, insertstep, nalleles, rnatype)
    genetrap(root, r, label, gene, transitions, G, R, S, insertstep, nalleles, rnatype, fittedparam, tempfish, priorcv, propcv, onstates)
end

function genetrap(root, r, label::String, gene::String, transitions::Tuple, G::Int, R::Int, S::Int, insertstep::Int, nalleles::Int=2, rnatype::String="", fittedparam=collect(1:num_rates(transitions, R,S,insertstep)-1), tempfish=1.0, priorcv=10.0, propcv=0.01, onstates=[])
    data = tempfish < 0 ? data_genetrap_FISH(root, label, gene) : data_genetrap(root, label, gene, tempfish)
    model = model_genetrap(gene, r, transitions, G, R, S, insertstep, fittedparam, nalleles, data.nRNA + 2, priorcv, propcv, onstates, rnatype)
    return data, model
end

function data_genetrap(root, label, gene, tempfish=1.0)
    LC = readLCPDF_genetrap(root, gene)
    if tempfish == 0
        counts = Int(div(sum(LC[:, 2] + LC[:, 3]), 2))
        println(counts)
        # set FISH counts to mean of live cell counts
        histFISH = readFISH_genetrap(root, gene, counts)
    else
        histFISH = readFISH_genetrap(root, gene, tempfish)
    end
    RNALiveCellData(label, gene, length(histFISH), histFISH, LC[:, 1], LC[:, 3], LC[:, 2])
end

function data_genetrap_FISH(root, label, gene)
    histFISH = readFISH_genetrap(root, gene, 1.0)
    RNAData(label, gene, length(histFISH), histFISH)
end

"""
model_genetrap

load model structure
"""

function model_genetrap(gene::String, r, transitions, G::Int, R::Int, S::Int, insertstep::Int, fittedparam::Vector, nalleles::Int, nhist::Int, priorcv=10.0, propcv=0.01, onstates=Int[], rnatype="", method=1)
    if S > 0
        S = R
    end
    if isempty(onstates)
        onstates = on_states(G, R, S, insertstep)
    end
    rm = fill(0.1, num_rates(transitions, R, S, insertstep))
    rcv = [fill(priorcv, length(rm) - 1); 0.1]
    if gene ∈ genes_gt()
        rm[end] = log(2.0) / (60 .* halflife_gt()[gene])
    end
    if r == 0.0
        r = rm
    end
    println(r)
    d = distribution_array(log.(rm[fittedparam]), sigmalognormal(rcv[fittedparam]), Normal)
    components = make_components_MTAI(transitions, G, R, S, insertstep, on_states(G, R, S, insertstep), nhist, r[num_rates(transitions, R, S, insertstep)])
    return GRSMmodel{typeof(r),typeof(d),typeof(propcv),typeof(fittedparam),typeof(method),typeof(components),typeof(onstates)}(G, R, S, insertstep, nalleles, rnatype, r, d, propcv, fittedparam, method, transitions, components, onstates)
end

# """
# setpriorrate(G,R,gene,Gprior::Tuple,Sprior::Tuple,Rcv::Float64)
# Set prior distribution for mean and cv of rates
# """
# function prior_rates_genetrap(G, R, gene, Gprior::Tuple, Sprior::Tuple, Rcv::Float64)
#     n = G - 1
#     rm = [fill(Gprior[1], 2 * n); .1; fill(.1, R); fill(Sprior[1], R); log(2.0) / (60 .* halflife_gt()[gene])]
#     rcv = [fill(Gprior[2], 2 * n); Rcv; Rcv; fill(Sprior[2], 2 * R); 0.01]
#     return rm, rcv
# end


"""
readrates_genetrap(infolder::String,rtype::String,gene::String,label,G,R,nalleles,rnatype::String)

Read in initial rates from previous runs
"""

function readrates_genetrap(infolder::String, rtype::String, gene::String, label, G, R, S, insertstep, nalleles, rnatype::String)
    row = get_row()[rtype]
    if rnatype == "offeject" || rnatype == "on"
        rnatype = ""
    end
    readrates_genetrap(getratefile_genetrap(infolder, rtype, gene, label, G, R, S, insertstep, nalleles, rnatype), row)
end

function readrates_genetrap(infile::String, row::Int)
    if isfile(infile) && ~isempty(read(infile))
        println(infile, ", row: ", get_rtype()[row])
        return readrates(infile, row, true)
    else
        println("using default rates")
        return 0
    end
end
"""
readLCPDF_genetrap(root,gene)

Read in dwell time PDF
"""
function readLCPDF_genetrap(root, gene)
    infile = joinpath(root, "DwellTimePDF/$(gene)_PDF.csv")
    if isfile(infile)
        LC = readdlm(infile, ',')
        x = truncate_histogram(LC[:, 2], 0.999, 1000)
        LC[1:length(x), :]
    else
        println("No data for gene: ", gene)
    end
end

"""
readFISH_genetrap(root,gene,temp=1.,clone=true)

Read and assemble smFISH data
"""
function readFISH_genetrap(root::String, gene::String, temp::Float64=1.0, clone=true)
    counts = countsFISH_genetrap(root, gene, clone)
    readFISH_genetrap(root, gene, Int(div(counts, temp)), clone)
end

function readFISH_genetrap(root::String, gene::String, counts::Int, clone=true)
    fishfile = joinpath(root, "Fish_4_Genes/$(gene)_steady_state_mRNA.csv")
    col = clone ? 3 : 2
    # Get smFISH RNA histograms
    if isfile(fishfile)
        x = readdlm(fishfile, ',')[3:end, col]
        x = x[x.!=""]
        x = truncate_histogram(x, 0.99, 1000)
        x /= sum(x)
        x *= counts
        return x
    else
        println("No smFISH data for gene: ", gene)
    end
    nothing
end

"""
countsFISH_genetrap(root,gene,x,counts,clone=true)

read in smFISH counts
"""
function countsFISH_genetrap(root, gene, clone=true)
    countfile = joinpath(root, "counts/total_counts.csv")
    wttext = "WT_"
    clonetext = "_clone"
    if isfile(countfile)
        countsdata = readdlm(countfile, ',')[:, :]
        if clone
            counts = (countsdata[findall(uppercase("$(gene)$clonetext") .== uppercase.(countsdata[:, 1])), 2])[1]
        else
            counts = (countsdata[findall("$wttext$(gene)" .== countsdata[:, 1]), 2])[1]
        end
    else
        counts = 1000
    end
    return counts
end

"""
survival_fraction(nu,eta,nr)
Fraction of introns that are not spliced prior to ejection
"""
function survival_fraction(nu, eta, nr)
    pd = 1.0
    for i = 1:nr
        pd *= nu[i+1] / (nu[i+1] + eta[i])
    end
    return pd
end


# """
# get_gamma(r,n,nr)
# G state forward and backward transition rates
# for use in transition rate matrices of Master equation
# (different from gamma used in Gillespie algorithms)
# """
# function get_gamma(r, n)
#     gammaf = zeros(n + 2)
#     gammab = zeros(n + 2)
#     for i = 1:n
#         gammaf[i+1] = r[2*(i-1)+1]
#         gammab[i+1] = r[2*i]
#     end
#     return gammaf, gammab
# end
# """
# get_nu(r,n,nr)
# R step forward transition rates
# """
# function get_nu(r, n, nr)
#     r[2*n+1:2*n+nr+1]
# end
# """
# get_eta(r,n,nr)
# Intron ejection rates at each R step
# """
# function get_eta(r, n, nr)
#     eta = zeros(nr)
#     if length(r) > 2 * n + 2 * nr
#         eta[1] = r[2*n+1+nr+1]
#         for i = 2:nr
#             # eta[i] = eta[i-1] + r[2*n + 1 + nr + i]
#             eta[i] = r[2*n+1+nr+i]
#         end
#     end
#     return eta
# end


# """
# Parameter information for GR models
# """
# function num_rates(transitions::Tuple, R)
#     length(transitions) + 2 * R + 2
# end
# function num_rates(G::Int, R)
#     2 * G + 2 * R
# end
# function Grange(G::Int)
#     n = G - 1
#     1:2*n
# end
# function initiation(G::Int)
#     n = G - 1
#     2 * n + 1
# end
# function Rrange(G::Int, R)
#     n = G - 1
#     2*n+2:2*n+R+1
# end
# function Srange(G::Int, R)
#     n = G - 1
#     2*n+R+2:2*n+2*R+1
# end

# """
# priordistributionLogNormal_genetrap(r,cv,G,R)

# LogNormal Prior distribution
# """
# function priordistributionLogNormal_genetrap(r, cv, G, R)
#     sigma = sigmalognormal(cv)
#     d = []
#     j = 1
#     #  G rates
#     for i in Grange(G)
#         push!(d, Distributions.LogNormal(log(r[i]), sigma[j]))
#         j += 1
#     end
#     # initiation rate
#     i = initiation(G)
#     push!(d, Distributions.LogNormal(log(r[i]), sigma[j]))
#     j += 1
#     # Step rates are bounded by length of insertion site to end of gene, i.e sum of reciprocal rates is bounded
#     t = sum(1 ./ r[Rrange(G, R)])
#     push!(d, Distributions.LogNormal(-log(t), sigma[j]))
#     j += 1
#     # priors for splice rates
#     rs = 0.0
#     for i in Srange(G, R)
#         rs = r[i]
#         push!(d, Distributions.LogNormal(log(rs), sigma[j]))
#         j += 1
#     end
#     return d
# end

# """
# priordistributionGamma_genetrap(r,cv,G,R)

# Gamma Prior distribution
# """
# function priordistributionGamma_genetrap(rm, cv, G, R)
#     n = G - 1
#     zet = R
#     d = []
#     rsig = cv .^ 2
#     j = 1
#     # priors for G rates
#     for i in Grange(G)
#         theta = rsig[j] * rm[i]
#         alpha = 1 / rsig[j]
#         push!(d, Distributions.Gamma(alpha, theta))
#         j += 1
#     end
#     # prior for nu1 is upper bounded by rm, i.e. distance is bounded by distance to insertion site
#     theta = rsig[j] * rm[initiation(G)]
#     alpha = 1 / rsig[j]
#     push!(d, Distributions.Gamma(alpha, theta))
#     j += 1
#     # priors for nu rates are bounded by length of insertion site to end of gene, i.e sum of reciprocal rates is bounded
#     tm = sum(1 ./ rm[Rrange(G, R)])
#     # for i in 2*n + 2 : 2*n + zet + 1
#     #     tm += 1/rm[i]
#     # end
#     theta = rsig[j] / tm
#     alpha = 1 / rsig[j]
#     push!(d, Distributions.Gamma(alpha, theta))
#     j += 1
#     # priors for ejection rates
#     rs = 0
#     for i in Srange(G, R)
#         rs = rm[i]
#         theta = rsig[j] * rs
#         alpha = 1 / rsig[j]
#         push!(d, Distributions.Gamma(alpha, theta))
#         j += 1
#     end
#     return d
# end
