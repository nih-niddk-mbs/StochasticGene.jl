# This file is part of StochasticGene.jl
#
#genetrap.jl
#Fit GRSM models to live cell and smFISH data


"""
Gene information
"""
genes_gt() = ["CANX";"DNAJC5";"ERRFI1";"KPNB1";"MYH9";"Rab7a";"RHOA";"RPAP3";"Sec16A";"SLC2A1"]
genelength_gt() = Dict([("Sec16A", 42960);("SLC2A1", 33802);("ERRFI1", 14615);("RHOA", 52948);("KPNB1", 33730);("MYH9", 106741);("DNAJC5", 40930);("CANX", 32710);("Rab7a", 88663);("RPAP3", 44130)])
MS2end_gt() = Dict([("Sec16A", 5220);("SLC2A1", 26001);("ERRFI1", 5324);("RHOA", 51109);("KPNB1", 24000);("MYH9", 71998);("DNAJC5", 14857);("CANX", 4861);("Rab7a", 83257);("RPAP3", 38610)])
halflife_gt() = Dict([("CANX", 50.),("DNAJC5", 5.),("ERRFI1", 1.35),("KPNB1", 9.),("MYH9", 10.),("Rab7a", 50.),("RHOA", 50.),("RPAP3", 7.5),("Sec16A", 8.),("SLC2A1", 5.)])


"""
    fit_genetrap()

"""

function fit_genetrap(nchains,maxtime,gene::String,transitions,G::Int,R::Int,insertstep::Int,infolder::String,folder::String,samplesteps::Int=1000;nalleles::Int=2,label="gt",rnatype="",fittedparam=collect(1:num_rates(transitions,R)-1),warmupsteps=0,annealsteps=0,temp=1.,tempanneal=100.,tempfish=1.,root::String=".",burst=false)
    println(now())
    data,model = genetrap(root,gene,transitions,G,R,insertstep,2,rnatype,fittedparam,infolder,folder,label,"median",tempfish)
    println("size of histogram: ",data.nRNA)
    options = MHOptions(samplesteps,warmupsteps,annealsteps,maxtime,temp,tempanneal)
    print_ll(data,model)
    fit,stats,measures = run_mh(data,model,options,nchains)
    if burst
        bs = burstsize(fit,model)
    else
        bs = 0
    end
    # finalize(data,model,fit,stats,measures,1.,folder,0,bs,root)
    finalize(data,model,fit,stats,measures,temp,folder,0,burst,false,root)
    println(now())
    return data, model_genetrap(get_rates(fit.parml,model),gene,G,R,insertstep,nalleles,fittedparam,rnatype,1,transitions,data)
end

"""
genetrap()

Load data, model and option structures for metropolis-hastings fit
root is the folder containing data and results

FISH counts is divided by tempfish to adjust relative weights of data
set tempfish = 0 to equalize FISH and live cell counts
"""
function genetrap(root,gene::String,transitions::Tuple,G::Int,R::Int,insertstep::Int,nalleles,rnatype::String,fittedparam::Vector,infolder::String,resultfolder::String,label::String,rtype::String,tempfish)
    # r = readrates_genetrap(root,infolder,rtype,gene,"$G$R",rnatype)
    r = readrates_genetrap(infolder,rtype,gene,label,G,R,nalleles,rnatype)
    println(r)
    genetrap(root,r,label,gene,G,R,transitions,insertstep,nalleles,rnatype,fittedparam,tempfish)
end

function genetrap(root,r,label::String,gene::String,G::Int,R::Int,transitions::Tuple,insertstep::Int,nalleles::Int=2,rnatype::String="",fittedparam=collect(1:num_rates(transitions,R)-1),tempfish=1.,method=1)
    # data = iszero(R) || iszero(tempfish) ? data_genetrap_FISH(root,label,gene) : data_genetrap(root,label,gene,tempfish)
    data = iszero(tempfish) ? data_genetrap_FISH(root,label,gene) : data_genetrap(root,label,gene,tempfish)
    model = model_genetrap(r,gene,G,R,insertstep,nalleles,fittedparam,rnatype,method,transitions,data.nRNA+2)
    return data,model
end

function data_genetrap(root,label,gene,tempfish=1.)
    LC = readLCPDF_genetrap(root,gene)
    if tempfish == 0
        counts = Int(div(sum(LC[:,2]+LC[:,3]),2))
        println(counts)
        # set FISH counts to mean of live cell counts
        histFISH = readFISH_genetrap(root,gene,counts)
    else
        histFISH = readFISH_genetrap(root,gene,tempfish)
    end
    RNALiveCellData(label,gene,length(histFISH),histFISH,LC[:,1],LC[:,3],LC[:,2])
end

function data_genetrap_FISH(root,label,gene)
    histFISH = readFISH_genetrap(root,gene,1.)
    RNAData(label,gene,length(histFISH),histFISH)
end

"""
model_genetrap

load model structure
"""

function model_genetrap(data,gene::String,transitions::Tuple,G::Int,R::Int,insertstep::Int,nalleles,rnatype::String,fittedparam::Vector,fixedeffects,infolder::String,label::String,rtype::String,root::String)
    # r = readrates_genetrap(root,infolder,rtype,gene,"$G$R",rnatype)
    r = readrates_genetrap(joinpath(root,infolder),rtype,gene,label,G,R,nalleles,rnatype)
    println(r)
    # genetrap(root,r,label,gene,G,R,transitions,nalleles,rnatype,fittedparam,tempfish)
    model_genetrap(r,gene,G,R,insertstep,nalleles,fittedparam,rnatype,method,transitions,data.nRNA+2)
end

function model_genetrap(r,gene::String,G::Int,R::Int,insertstep::Int,nalleles::Int,fittedparam,rnatype::String,method,transitions,nhist)
    ntransitions = length(transitions)
    S = R
    if R == 0
        decayrate = get_decay(halflife_gt()[gene])
        if r == 0
            r = setrate(G,1,decayrate,.1)[1]
        end
        d = prior_rna(r,G,1,fittedparam,decayrate,1.)
        components = make_components(transitions,G,r,nhist,Indices(collect(1:ntransitions),[ntransitions+1],Int[],ntransitions + 2))
        return GMmodel{typeof(r),typeof(d),Float64,typeof(fittedparam),typeof(method),typeof(components)}(G,nalleles,r,d,0.05,fittedparam,method,transitions,components,onstates)
    else
        propcv=[.02*ones(2*(G-1));.02;.02*ones(R);.02*ones(R)]
        Gprior=(.01,10)
        Sprior=(.1,10)
        Rcv=10.
        if gene ∈ genes_gt()
            rm,rcv = prior_rates_genetrap(G,R,gene,Gprior,Sprior,Rcv)
        else
            rm =[fill(0.01, length(rtarget) - 1);r[end]]
            rcv = fill(cv, length(r))
        end
        if r == 0
            r = rm
            println(r)
        end
        d = distribution_array(log.(rm[fittedparam]),sigmalognormal(rcv[fittedparam]),Normal)
        components = make_components_MTAI(transitions, G, R, S, insertstep, on_states(G, R, S, insertstep), nhist, r[num_rates(transitions,R,S,insertstep)])
        return GRSMmodel{typeof(r),typeof(d),typeof(propcv),typeof(fittedparam),typeof(method),typeof(components),Int}(G,R,1,insertstep,nalleles,rnatype,r,d,propcv,fittedparam,method,transitions,components,0)
    end
end

"""
setpriorrate(G,R,gene,Gprior::Tuple,Sprior::Tuple,Rcv::Float64)
Set prior distribution for mean and cv of rates
"""
function prior_rates_genetrap(G,R,gene,Gprior::Tuple,Sprior::Tuple,Rcv::Float64)
    n = G-1
    rm = [fill(Gprior[1],2*n);2000/(genelength_gt()[gene]-MS2end_gt()[gene]);fill(2000/MS2end_gt()[gene]*R,R);fill(Sprior[1], R);log(2.)/(60 .* halflife_gt()[gene])]
    rcv = [fill(Gprior[2],2*n);Rcv;Rcv;fill(Sprior[2], 2*R);.01]
    return rm,rcv
end

"""
priordistributionLogNormal_genetrap(r,cv,G,R)

LogNormal Prior distribution
"""
function priordistributionLogNormal_genetrap(r,cv,G,R)
    sigma = sigmalognormal(cv)
    d = []
    j = 1
    #  G rates
    for i in Grange(G)
        push!(d,Distributions.LogNormal(log(r[i]),sigma[j]))
        j += 1
    end
    # initiation rate
    i = initiation(G)
    push!(d,Distributions.LogNormal(log(r[i]),sigma[j]))
    j += 1
    # Step rates are bounded by length of insertion site to end of gene, i.e sum of reciprocal rates is bounded
    t = sum(1 ./ r[Rrange(G,R)])
    push!(d,Distributions.LogNormal(-log(t),sigma[j]))
    j += 1
    # priors for splice rates
    rs = 0.
    for i in Srange(G,R)
        rs = r[i]
        push!(d,Distributions.LogNormal(log(rs),sigma[j]))
        j += 1
    end
    return d
end

"""
priordistributionGamma_genetrap(r,cv,G,R)

Gamma Prior distribution
"""
function priordistributionGamma_genetrap(rm,cv,G,R)
    n = G-1
    zet = R
    d = []
    rsig = cv .^ 2
    j = 1
    # priors for G rates
    for i in Grange(G)
        theta = rsig[j]*rm[i]
        alpha = 1/rsig[j]
        push!(d,Distributions.Gamma(alpha,theta))
        j += 1
    end
    # prior for nu1 is upper bounded by rm, i.e. distance is bounded by distance to insertion site
    theta = rsig[j]*rm[initiation(G)]
    alpha = 1/rsig[j]
    push!(d,Distributions.Gamma(alpha,theta))
    j += 1
    # priors for nu rates are bounded by length of insertion site to end of gene, i.e sum of reciprocal rates is bounded
    tm = sum(1 ./ rm[Rrange(G,R)])
    # for i in 2*n + 2 : 2*n + zet + 1
    #     tm += 1/rm[i]
    # end
    theta = rsig[j]/tm
    alpha = 1/rsig[j]
    push!(d,Distributions.Gamma(alpha,theta))
    j += 1
    # priors for ejection rates
    rs = 0
    for i in Srange(G,R)
        rs = rm[i]
        theta = rsig[j]*rs
        alpha = 1/rsig[j]
        push!(d,Distributions.Gamma(alpha,theta))
        j += 1
    end
    return d
end


"""
readLCPDF_genetrap(root,gene)

Read in dwell time PDF
"""
function readLCPDF_genetrap(root,gene)
    infile = joinpath(root,"DwellTimePDF/$(gene)_PDF.csv")
    if isfile(infile)
        LC = readdlm(infile,',')
        x = truncate_histogram(LC[:,2],.999,1000)
        LC[1:length(x),:]
    else
        println("No data for gene: ", gene)
    end
end

"""
readFISH_genetrap(root,gene,temp=1.,clone=true)

Read and assemble smFISH data
"""
function readFISH_genetrap(root::String,gene::String,temp::Float64=1.,clone=true)
    counts = countsFISH_genetrap(root,gene,clone)
    readFISH_genetrap(root,gene,Int(div(counts,temp)),clone)
end

function readFISH_genetrap(root::String,gene::String,counts::Int,clone=true)
    fishfile=joinpath(root,"Fish_4_Genes/$(gene)_steady_state_mRNA.csv")
    col = clone ? 3 : 2
    # Get smFISH RNA histograms
    if isfile(fishfile)
        x = readdlm(fishfile,',')[3:end,col]
        x = x[x .!= ""]
        x = truncate_histogram(x,.99,1000)
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
function countsFISH_genetrap(root,gene,clone=true)
    countfile = joinpath(root,"counts/total_counts.csv")
    wttext = "WT_"
    clonetext = "_clone"
    if isfile(countfile)
        countsdata = readdlm(countfile,',')[:,:]
        if clone
            counts=(countsdata[findall(uppercase("$(gene)$clonetext") .== uppercase.(countsdata[:,1])),2])[1]
        else
            counts=(countsdata[findall("$wttext$(gene)" .== countsdata[:,1]),2])[1]
        end
    else
        counts = 1000
    end
    return counts
end

"""
get_gamma(r,n,nr)
G state forward and backward transition rates
for use in transition rate matrices of Master equation
(different from gamma used in Gillespie algorithms)
"""
function get_gamma(r,n)
    gammaf = zeros(n+2)
    gammab = zeros(n+2)
    for i = 1:n
        gammaf[i+1] = r[2*(i-1)+1]
        gammab[i+1] = r[2*i]
    end
    return gammaf, gammab
end
"""
get_nu(r,n,nr)
R step forward transition rates
"""
function get_nu(r,n,nr)
    r[2*n+1 : 2*n+nr+1]
end
"""
get_eta(r,n,nr)
Intron ejection rates at each R step
"""
function get_eta(r,n,nr)
    eta = zeros(nr)
    if length(r) > 2*n + 2*nr
        eta[1] = r[2*n + 1 + nr + 1]
        for i = 2:nr
            # eta[i] = eta[i-1] + r[2*n + 1 + nr + i]
            eta[i] = r[2*n + 1 + nr + i]
        end
    end
    return eta
end
"""
survival_fraction(nu,eta,nr)
Fraction of introns that are not spliced prior to ejection
"""
function survival_fraction(nu,eta,nr)
    pd = 1.
    for i = 1:nr
        pd *= nu[i+1]/(nu[i+1]+eta[i])
    end
    return pd
end


"""
Parameter information for GR models
"""
function num_rates(transitions::Tuple,R)
    length(transitions) + 2*R + 2
end
function num_rates(G::Int,R)
    2*G + 2*R
end
function Grange(G::Int)
    n = G - 1
    1 : 2*n
end
function initiation(G::Int)
    n = G - 1
    2*n + 1
end
function Rrange(G::Int,R)
    n = G - 1
    2*n + 2 : 2*n + R + 1
end
function Srange(G::Int,R)
    n = G - 1
    2*n + R + 2 : 2*n + 2*R + 1
end
