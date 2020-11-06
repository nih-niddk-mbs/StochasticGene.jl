
"""
genetrap.jl

Fit GRSM models to live cell and smFISH data
"""

"""
Gene information
"""
const genes = ["CANX";"DNAJC5";"ERRFI1";"KPNB1";"MYH9";"Rab7a";"RHOA";"RPAP3";"Sec16A";"SLC2A1"]
const genelength = Dict([("Sec16A", 42960);("SLC2A1", 33802);("ERRFI1", 14615);("RHOA", 52948);("KPNB1", 33730);("MYH9", 106741);("DNAJC5", 40930);("CANX", 32710);("Rab7a", 88663);("RPAP3", 44130)])
const MS2end = Dict([("Sec16A", 5220);("SLC2A1", 26001);("ERRFI1", 5324);("RHOA", 51109);("KPNB1", 24000);("MYH9", 71998);("DNAJC5", 14857);("CANX", 4861);("Rab7a", 83257);("RPAP3", 38610)])
const halflife = Dict([("CANX", 50.),("DNAJC5", 5.),("ERRFI1", 1.35),("KPNB1", 9.),("MYH9", 10.),("Rab7a", 50.),("RHOA", 50.),("RPAP3", 7.5),("Sec16A", 8.),("SLC2A1", 5.)])
const models = [21;22;23;34;35;31;32;33;34;35]
const datapath = "/Users/carsonc/Dropbox/Larson/GeneTrap_analysis"
const resultpath = "/Users/carsonc/Dropbox/Larson/GeneTrap_analysis/Results"

function write_results(fit,model::GRSMmodel)
    println(fit[1].llml)
end

function genetrap(infolder::String,rinchar::String,gene::String,G::Int,R::Int,nalleles::Int,type::String,maxtime::Float64,samplesteps::Int,annealsteps=0,warmupsteps=0,temp=1)
    data = data_genetrap(gene)
    model = model_genetrap(infolder,rinchar,gene,G,R,nalleles,type)
    options = MHOptions(samplesteps,annealsteps,warmupsteps,maxtime,temp)
    return data,model,options
end

function data_genetrap(gene)
    LC = readLCPDF_genetrap(gene)
    histFISH = readFISH_genetrap(gene)
    RNALiveCellData(gene,LC[:,1],LC[:,2],LC[:,3],length(histFISH),histFISH)
end

function model_genetrap(infolder::String,rinchar::String,gene::String,G::Int,R::Int,nalleles::Int,type::String,propcv=[.05*ones(2*(G-1));1e-4;1e-4*ones(R);.05*ones(R)],fittedparam=collect(1:num_rates(G,R)-1),method=0,Gprior=(.01,.5),Sprior=(.1,10.),Rcv=1e-2)
    rm,rcv = setpriorrate(G,R,gene,Gprior,Sprior,Rcv)
    d = priordistributionLogNormal(rm,rcv,G,R)
    r = read_rates(infolder,rinchar,gene,"$G$R",type)[:,1]
    if r == 0
        r = rm
    end
     GRSMmodel{typeof(r),typeof(d),typeof(propcv),typeof(fittedparam),typeof(method)}(G,R,nalleles,type,r,d,propcv,fittedparam,method)
    # GRSMmodel(G,R,nalleles,type,r,d,prop,fittedparam,method)
end

"""
likelihoodfn(r,data::RNALiveCellData,model::GRSMmodel)
called by loglikelihood in metropolis_hastings.jl
"""
# function likelihoodfn(param,data::RNALiveCellData,model::GRSMmodel)
function likelihoodfn(param::Vector,data::RNALiveCellData,model::GRSMmodel,method)
    pdftuple = likelihoodtuple(param,data,model)
    pdfvector = Array{Float64,1}(undef,0)
    for p in pdftuple
        vcat(pdfvector,p)
    end
    return pdfvector
end

# function likelihoodfn(param,data::RNALiveCellData,model::GRSMmodel)
function likelihoodfn(param::Vector,data::RNALiveCellData,model::GRSMmodel)
    modelOFF, modelON, histF = likelihoodtuple(param,data,model)
    return [modelOFF;modelON;histF]
end

function likelihoodtuple(param,data::RNALiveCellData,model::GRSMmodel)
    r = get_rates_genetrap(param,model)
    if model.method == 0
        modelOFF, modelON = offonPDF(data.bins,r,model.G-1,model.R)
        if model.type == "off"
            histF = steady_state_offpath(r,model.G-1,model.R,data.nRNA,model.nalleles)
        else
            histF = steady_state(r,model.G-1,model.R,data.nRNA,model.nalleles)
        end
    else
        if model.type == "off"
            modelOFF,modelON,histF = telegraphoff(data.bins,data.nRNA,r,model.G-1,model.R,model.nalleles)
        else
            modelOFF,modelON,histF = telegraph(data.bins,data.nRNA,r,model.G-1,model.R,model.nalleles)
        end
    end
    return modelOFF, modelON, histF
end

function get_rates(params,model::GRSMmodel)
    r = copy(model.rates)
    r[model.fittedparams] = params
    return r
end


"""
logprior(x,model::GRSModel)
Compute log prior using distributions in Model.rateprior
called by mhstep() in metropolis_hastings.jl
"""
function logprior(param,model::GRSMmodel)
    r = get_rates(param,model)
    d = model.rateprior
    G = model.G
    R = model.R
    p=0
    j = 1
    #  G rates
    for i in Grange(G)
        p -= logpdf(d[j],r[i])
        j += 1
    end
    # initiation rate
    i = initiation(G)
    p -= logpdf(d[j],r[i])
    j += 1
    # sum of R Steps rates are bounded by length of insertion site to end of gene, i.e sum of reciprocal rates is bounded
    t = sum(1 ./ r[Rrange(G,R)])
    p -= logpdf(d[j],1/t)
    j += 1
    # priors for splice rates
    rs = 0
    for i in Srange(G,R)
        rs += r[i]
        p -= logpdf(d[j],rs)
        j += 1
    end
    return p
end

"""
setpriorrate(G,R,gene,Gprior::Tuple,Sprior::Tuple,Rcv::Float64)
Set prior distribution for mean and cv of rates
"""
function setpriorrate(G,R,gene,Gprior::Tuple,Sprior::Tuple,Rcv::Float64)
    n = G-1
    rm = [fill(Gprior[1],2*n);2000/(genelength[gene]-MS2end[gene]);fill(2000/MS2end[gene]*R,R);fill(Sprior[1], R);log(2.)/(60 .* halflife[gene])]
    rcv = [fill(Gprior[2],2*n);Rcv;Rcv;fill(Sprior[2], R);Rcv]
    return rm,rcv
end

"""
priordistributionLogNormal(r,cv,G,R)
LogNormal Prior distribution
"""
function priordistributionLogNormal(r,cv,G,R)
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
        rs += r[i]
        push!(d,Distributions.LogNormal(log(rs),sigma[j]))
        j += 1
    end
    return d
end

function priordistributionGamma(rm,cv,G,R)
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
        rs += rm[i]
        theta = rsig[j]*rs
        alpha = 1/rsig[j]
        push!(d,Distributions.Gamma(alpha,theta))
        j += 1
    end
    return d
end

function datahistogram(data::RNALiveCellData)
    return [data.OFF;data.ON;data.histRNA]
end

function likelihoodtuple(r,data::RNALiveCellData,model,method)
    telegraph(data.bins,model.nRNA,n,nr,r,model.nalleles)
end

function likelihoodtuples(bins,nRNA,r,n,nr,nalleles,type)
    if type == "off"
        modelOFFt, modelONt, histFt = telegraphoff(bins,nRNA,r,n,nr,nalleles)
        histF = steady_state_offpath(r,n,nr,nRNA,nalleles)
    else
        modelOFFt, modelONt, histFt = telegraph(bins,nRNA,r,n,nr,nalleles)
        histF = steady_state(r,n,nr,nRNA,nalleles)
    end
    modelOFF, modelON = offonPDF(bins,r,n,nr)
    return modelOFF,modelON,histF,modelOFFt,modelONt,histFt
end

"""
Parameter information
"""
function num_rates(G,R)
    n = G - 1
    2*n + 2*R + 2
end
function Grange(G)
    n = G - 1
    1 : 2*n
end
function initiation(G)
    n = G - 1
    2*n + 1
end
function Rrange(G,R)
    n = G - 1
    2*n + 2 : 2*n + R + 1
end
function Srange(G,R)
    n = G - 1
    2*n + R + 2 : 2*n + 2*R + 1
end
"""
Read in dwell time PDF
"""
function readLCPDF_genetrap(gene)
    infile = "$datapath/DwellTimePDF/$(gene)_PDF.csv"
    if isfile(infile)
        LC = readdlm(infile,',')
        x = truncate_histogram(LC[:,2],.999,1000)
        LC[1:length(x),:]
    else
        println("No data for gene: ", gene)
    end
end

"""
Read and assemble smFISH data
"""
function readFISH_genetrap(gene,Clone=true)
    countfile="$datapath/counts/total_counts.csv"
    fishfile="$datapath/Fish_4_Genes/$(gene)_steady_state_mRNA.csv"
    wttext = "WT_"
    clonetext = "_clone"
    if Clone
        col = 3
    else
        col = 2
    end
    # Get total RNA counts
    if isfile(countfile)
        countsdata = readdlm(countfile,',')[:,:]
        if Clone
            counts=(countsdata[findall(uppercase("$(gene)$clonetext") .== uppercase.(countsdata[:,1])),2])[1]
        else
            counts=(countsdata[findall("$wttext$(gene)" .== countsdata[:,1]),2])[1]
        end
    else
        counts = 2000
    end
    # Get smFISH RNA histograms
    if isfile(fishfile)
        x = readdlm(fishfile,',')[3:end,col]
        x = x[x .!= ""]
        x /= sum(x)
        x *= counts
        return truncate_histogram(x,.99,1000)
    else
        println("No smFISH data for gene: ", gene)
    end
    nothing
end

"""
read_rates(infolder::String,rinchar::String,gene::String,model::String,type::String)
Read in initial rates from previous runs
"""
function read_rates_genetrap(infolder::String,rinchar::String,gene::String,model::String,type::String)
    if rinchar == "rml" || rinchar == "rlast"
        rskip = 1
        rstart = 1
    elseif rinchar == "rmean"
        rskip = 2
        rstart = 1
    elseif rinchar == "rquant"
        rskip = 3
        rstart = 2
    end
    infile = "$resultpath/$infolder/$(gene)/$rinchar$model$type$txtstr"
    if isfile(infile) && ~isempty(read(infile))
        r = readdlm(infile)
        r = r[:,rstart:rskip:end]
        println(gene," ",model," ",type)
        return r
    else
        println(gene," ",model," ",type," no prior")
        return 0
    end
    nothing
end
