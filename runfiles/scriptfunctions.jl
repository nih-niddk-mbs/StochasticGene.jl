
using Dates
using DelimitedFiles


function fit_rna(nchains,gene::String,cond::String,G::Int,maxtime::Float64,infolder::String,resultfolder,datafolder,inlabel,label,nsets,runcycle::Bool=false,sample::Int=40000,warmup=20000,anneal=100000,temp=1.,tempanneal=100.,root = "/home/carsonc/scrna/")

    params = (sample,warmup,anneal,temp,tempanneal)

    if nsets == 1
        if G == 3
            fittedparam = [1,2,3,4,5,7]
        else
            fittedparam = [1,2,3,5]
        end
    else
        if G == 3
            fittedparam = 7 .+ [1,2,3,4,5,7]
        else
            fittedparam = 5 .+ [1,2,3,5]
        end
    end

    @time fit_rna(nchains,gene,cond,G,maxtime,fittedparam,infolder,resultfolder,datafolder,inlabel,label,nsets,runcycle,params,root)

end

function fit_rna(nchains::Int,gene::String,cond::String,G::Int,maxtime::Float64,fittedparam::Vector,infolder::String,resultfolder,datafolder,inlabel,label,nsets::Int,runcycle::Bool,params::Tuple,root)
    if cond == "null"
        cond = ""
    end
    if nsets == 1
        datafile = StochasticGene.scRNApath(gene,cond,datafolder,root)
        data = StochasticGene.data_rna(datafile,label,gene,false)
    else
        sets =["T0","T30","T120"]
        datafile =[]
        fittedparam .+= 2*G
        for set in sets
            folder = joinpath(datafolder,set)
            datafile = vcat(datafile,StochasticGene.scRNApath(gene,cond,folder,root))
        end
        println(datafile)
        # datafile = StochasticGene.scRNApath(gene,cond,datafolder,root)
        data = StochasticGene.data_rna(datafile,label,[0.,30.,120.],gene,false)
    end
    fit_rna(nchains,gene,G,data,maxtime,nsets,fittedparam,infolder,resultfolder,datafolder,inlabel,runcycle,params,root)
end

function fit_rna(nchains::Int,gene::String,G::Int,data::StochasticGene.HistogramData,maxtime::Float64,nsets,fittedparam,infolder,resultfolder,datafolder,inlabel,runcycle::Bool,params::Tuple,root)

    println(now())
    println(gene," ",G)
    println("in: ", infolder," out: ",resultfolder)
    println(datafolder)
    println(maxtime)

    # model parameters
    decayrate = decay(root,gene)
    println(decayrate)
    if decayrate < 0
        throw("error")
    end
    nalleles = alleles(root,gene)
    yieldprior = .1

    # Metropolis-Hastings parameters
    # maxtime = 3600. * 1/60
    samplesteps = params[1]
    warmupsteps = params[2]
    annealsteps = params[3]
    temp = params[4]
    tempanneal = params[5]

    # samplesteps = 200000
    # warmupsteps = 0
    # annealsteps = 0
    # temp = 200.
    # tempanneal = 200.

    r = getr(gene,G,nalleles,decayrate,fittedparam,inlabel,infolder,nsets,root)

    if runcycle
        r = cycle(nchains,data,r,G,nalleles,nsets,.02,fittedparam,decayrate,yieldprior,maxtime/2,temp,tempanneal)
    end

    cv = getcv(gene,G,nalleles,fittedparam,inlabel,infolder,root)
    options = StochasticGene.MHOptions(samplesteps,annealsteps,warmupsteps,maxtime/2,temp,tempanneal)
    model = StochasticGene.model_rna(r,G,nalleles,nsets,cv,fittedparam,decayrate,yieldprior,0)

    param,_ = StochasticGene.initial_proposal(model)
    ll,_ = StochasticGene.loglikelihood(param,data,model)
    println("initial ll: ",ll)

    fit,stats,waic = StochasticGene.run_mh(data,model,options,nchains);

    writefile = joinpath(root,resultfolder)
    StochasticGene.writeall(writefile,fit,stats,waic,data,temp,model)

    println("final ll: ",fit.llml)
    println(fit.accept," ",fit.total)
    println("Deviance: ",StochasticGene.deviance(fit,data,model))
    println(now())
end


function cycle(nchains,data,r,G,nalleles,nsets,cv,fittedparam,decayrate,yieldprior,maxtime,temp,tempanneal)
    options = StochasticGene.MHOptions(100,0,0,maxtime/10,temp,tempanneal)
    t0 = time()
    while (time() - t0 < maxtime)
        for fp in fittedparam
            model = StochasticGene.model_rna(r,G,nalleles,nsets,cv,[fp],decayrate,yieldprior,0)
            fit,stats,waic = StochasticGene.run_mh(data,model,options,nchains);
            r = StochasticGene.get_rates(fit.parml,model)
        end
    end
    return r
end

function getr(gene,G,nalleles,decayrate,fittedparam,inlabel,infolder,nsets,root)



end

function  getr(gene,G,nalleles,decayrate,fittedparam,inlabel,infolder,nsets,root)
    ratefile = StochasticGene.path_Gmodel("rates",gene,G,nalleles,inlabel,infolder,root)

    if isfile(ratefile)
        r = StochasticGene.readrates(ratefile,1)
        r[1:2*G] *= decayrate/r[2*G]
        r[2*G+1] = clamp(r[2*G+1],eps(Float64),1-eps(Float64))
    else
        if G == 2
            r = [0.015,0.015,0.5,.01,1.]*decayrate/.01
            r[end] = .1
        else
            r = [0.015,.2,.2,0.015,1.5,.01,1.]*decayrate/.01
            r[end] = .1
        end
    end
    if nsets == 2
        r = [r[1:2*G];r[1:2*G];r[2*G+1]]
    end
    println(r)
    return r

end

function getcv(gene,G,nalleles,fittedparam,inlabel,infolder,root)
    paramfile = StochasticGene.path_Gmodel("param_stats",gene,G,nalleles,inlabel,infolder,root)

    if isfile(paramfile)
        cv = StochasticGene.read_covlogparam(paramfile)
        cv = float.(cv)
        if ~StochasticGene.isposdef(cv) || size(cv)[1] != length(fittedparam)
            cv = .02
        end
    else
        cv = .02
    end
    println(cv)
    return cv
end

function decay(root::String,gene,file="data/HCT116_all_cells_histograms_and_half_lives_March_2021/Starved_Rep7_half_lives.csv",col=2)
    path = joinpath(root,file)
    if isfile(path)
        in = readdlm(path,',')
        a = in[findfirst(in[:,1] .== gene),col]
        println(a)
        return decay(a)
    end
    println(gene," has no decay time")
    return -1.
end

decay(a::Float64) = log(2)/a/60.

function decay(a,gene)
    if typeof(a) <: Number
        return decay(a)
    else
        println(gene," has no decay time")
        return -1.
    end
end

function alleles(root,gene,file="data/HCT116_alleles_number.txt")
    in = readdlm(joinpath(root,"data/HCT116_alleles_number.txt"))
    in[findfirst(in[:,1] .== gene),3]
end

function fix_measures(resultfolder,measurefile,ratefile::String,cond,n,datafolder,root)
    front,back = split(measurefile,".")
    fix = front * "fix" * "." * back
    measurepath = joinpath(resultfolder,measurefile)
    ratepath = joinpath(resultfolder,ratefile)
    rates,_ = readdlm(ratepath,',',header=true)
    println(length(rates[:,1]))
    measures,head = readdlm(measurepath,',',header=true)
    f = open(joinpath(resultfolder,fix),"w")
    writedlm(f,[head "Deviance fixed" "LogML fixed" "AIC fixed"],',')
    for i in eachindex(rates[:,1])
        h,hd = histograms(rates[i,:],cond,n,datafolder,root)
        d = StochasticGene.deviance(h,hd)
        ll = StochasticGene.crossentropy(h,hd)
        a = 4*(n+1) + 2*ll
        writedlm(f,[measures[i:i,:] d ll a],',')
    end
    close(f)
end

function compute_deviance(outfile,ratefile::String,cond,n,datafolder,root)
    f = open(outfile,"w")
    rates,head = readdlm(ratefile,',',header=true)
    for r in eachrow(rates)
        d=deviance(r,cond,n,datafolder,root)
        writedlm(f,[gene d],',')
    end
    close(f)
end

function write_histograms(outfolder,ratefile::String,cond,n,datafolder,root)
    rates,head = readdlm(ratefile,',',header=true)
    for r in eachrow(rates)
        h,hd = histograms(r,cond,n,datafolder,root)
        f = open(joinpath(outfolder,r[1] * ".txt"),"w")
        writedlm(f,h)
        close(f)
    end
end

function write_burst_stats(outfile,infile::String,n,root)
    f = open(outfile,"w")
    contents,head = readdlm(infile,',',header=true)
    writedlm(f,["Gene" "Mean OFF Period (min)" "Burst Size"],',')
    for r in eachrow(contents)
        gene = String(r[1])
        off = meanofftime(r[2:2*n+3],n,1)
        # h,hd = histograms(r,cond,n,datafolder,root)
        writedlm(f,[gene off r[2*n+1]/r[2*n]],',')
    end
    close(f)
end

meanofftime(gene::String,infile,n,method,root) = sum(1 .- offtime(gene,infile,n,method,root))

function meanofftime(r::Vector,n::Int,method::Int)
    if n == 1
        return 1/r[1]
    else
        return sum(1 .- offtime(r,n,method))
    end
end

function offtime(r::Vector,n::Int,method::Int)
    _,_,TI = StochasticGene.mat_G_DT(r,n)
    vals,_ = StochasticGene.eig_decompose(TI)
    minval = min(minimum(abs.(vals[vals.!=0])),.2)
    StochasticGene.offtimeCDF(collect(1.:5/minval),r,n,TI,method)
end

function offtime(gene::String,infile,n,method,root)
    contents,head = readdlm(infile,',',header=true)
    r = float.(contents[gene .== contents[:,1],2:end-1])[1,:]
    offtime(r,n,method)

end


# function offtime(gene,infile,n,method,root)
#     contents,head = readdlm(infile,',',header=true)
#     r = float.(contents[gene .== contents[:,1],2:end-1])[1,:]
#     println(r)
#     _,_,TI = StochasticGene.mat_G_DT(r,n)
#     vals,vects = StochasticGene.eig_decompose(TI)
#     minval = minimum(abs.(vals[vals.!=0]))
#     StochasticGene.offtimeCDF(collect(1.:5/minval),r,n,TI,method)
# end


function write_moments(outfile,infile::String,cond,datafolder,root)
    f = open(outfile,"w")
    contents,head = readdlm(infile,',',header=true)
    writedlm(f,["Gene" "Expression Mean" "Expression Variance"],',')
    for r in eachrow(contents)
        gene = String(r[1])
        datafile = StochasticGene.scRNApath(gene,cond,datafolder,root)
        # data = StochasticGene.data_rna(datafile,label,gene,false)
        len,h = StochasticGene.histograms_rna(datafile,gene,false)
        # h,hd = histograms(r,cond,n,datafolder,root)
        writedlm(f,[gene StochasticGene.mean_histogram(h) StochasticGene.var_histogram(h)],',')
    end
    close(f)
end

function histograms(r,cond,n,datafolder,root)
    gene = String(r[1])
    datafile = StochasticGene.scRNApath(gene,cond,datafolder,root)
    hd = StochasticGene.read_scrna(datafile)
    h = StochasticGene.steady_state(r[2:2*n+3],r[end],n,length(hd),alleles(root,gene))
    return h,hd
end

function deviance(r,cond,n,datafolder,root)
    h,hd = histograms(r,cond,n,datafolder,root)
    StochasticGene.deviance(h,hd)
end

function bestmodel(measures2,measures3)
    m2,head = readdlm(infile,',',header=true)
    m3,head = readdlm(infile,',',header=true)
end

function join_files(file1::String,file2::String,outfile::String,addlabel::Bool=true)
    contents1,head1 = readdlm(file1,',',header=true)   # model G=2
    contents2,head2 = readdlm(file2,',',header=true)   # model G=3
    f = open(outfile,"w")
    if addlabel
        header = vcat(String.(head1[2:end]) .* "_G2",String.(head2[2:end]) .* "_G3")
    else
        header = vcat(String.(head1[2:end]),String.(head2[2:end]))
    end
    header = reshape(permutedims(header),(1,length(head1)+length(head2)-2))
    header = hcat(head1[1],header)
    println(header)
    writedlm(f,header,',')
    for row in 1:size(contents1,1)
        if contents1[row,1] == contents2[row,1]
            contents = hcat(contents1[row:row,2:end],contents2[row:row,2:end])
            contents = reshape(permutedims(contents),(1,size(contents1,2)+size(contents2,2)-2))
            contents = hcat(contents1[row,1],contents)
            writedlm(f,contents,',')
        end
    end
    close(f)
end

function best_model(file::String)
    contents,head = readdlm(file,',',header=true)
    f = open(file,"w")
    head = hcat(head,"Winning Model")
    writedlm(f,head,',')
    for row in eachrow(contents)
        if abs(row[11] - row[4]) > row[5] + row[12]
            if row[11] < row[4]
                writedlm(f, hcat(permutedims(row),3),',')
            else
                writedlm(f, hcat(permutedims(row),2),',')
            end
        else
            if row[13] < row[6]
                writedlm(f, hcat(permutedims(row),3),',')
            else
                writedlm(f, hcat(permutedims(row),2),',')
            end
        end
    end
    close(f)
end

# function add_burst(filein,fileout,filemodel2,filemodel3)
#     contents,head = readdlm(filein,',',header=true)
#     burst2,head2 = readdlm(filemodel2,',',header=true)
#     burst3,head3 = readdlm(filemodel3,',',header=true)
#     f = open(fileout,"w")
#     head = hcat(head,["mean off period" "bust size"])
#     writedlm(f,head,',')
#     for row in eachrow(contents)
#         if Int(row[end]) == 2
#             writedlm(f, hcat(permutedims(row),permutedims(burst2[findfirst(burst2[:,1] .== row[1]),2:3])),',')
#         else
#             writedlm(f, hcat(permutedims(row),permutedims(burst3[findfirst(burst3[:,1] .== row[1]),2:3])),',')
#         end
#     end
#     close(f)
# end


function add_best_burst(filein,fileout,filemodel2,filemodel3)
    contents,head = readdlm(filein,',',header=true)
    burst2,head2 = readdlm(filemodel2,',',header=true)
    burst3,head3 = readdlm(filemodel3,',',header=true)
    f = open(fileout,"w")
    head = hcat(head,["mean off period" "bust size"])
    writedlm(f,head,',')
    for row in eachrow(contents)
        if Int(row[end]) == 2
            writedlm(f, hcat(permutedims(row),permutedims(burst2[findfirst(burst2[:,1] .== row[1]),2:3])),',')
        else
            writedlm(f, hcat(permutedims(row),permutedims(burst3[findfirst(burst3[:,1] .== row[1]),2:3])),',')
        end
    end
    close(f)
end

function add_best_occupancy(filein,fileout,filemodel2,filemodel3)
    contents,head = readdlm(filein,',',header=true)
    occupancy2,head2 = readdlm(filemodel2,',',header=true)
    occupancy3,head3 = readdlm(filemodel3,',',header=true)
    f = open(fileout,"w")
    head = hcat(head,["Off -2" "Off -1" "On State" ])
    writedlm(f,head,',')
    for row in eachrow(contents)
        if Int(row[end-2]) == 2
            writedlm(f, hcat(permutedims(row),hcat("NA",permutedims(occupancy2[findfirst(occupancy2[:,1] .== row[1]),2:end]))),',')
        else
            writedlm(f, hcat(permutedims(row),permutedims(occupancy3[findfirst(occupancy3[:,1] .== row[1]),2:end])),',')
        end
    end
    close(f)
end


# precompile(fit_rna,(String,String,Int,Float64,String,String,String,String,String,Int,String,Bool))
# fit_rna(nchains,gene::String,cond::String,G::Int,maxtime::Float64,infolder::String,resultfolder,datafolder,inlabel,label,nsets,runcycle::Bool=false,sample::Int=40000,warmup=20000,anneal=100000,temp=1.,tempanneal=100.,root = "/home/carsonc/scrna/")