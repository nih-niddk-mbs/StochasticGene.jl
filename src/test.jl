
using Distributions
using StatsBase
using Statistics
using DelimitedFiles
using Dates
using Distributed
using LinearAlgebra
using Plots
using SparseArrays
using DifferentialEquations
using LSODA
using DataFrames
using FFTW
using Downloads
using CSV
using MultivariateStats
using Optim

include("common.jl")
include("transition_rate_matrices_transpose.jl")
include("chemical_master_new.jl")
include("simulator.jl")
include("utilities.jl")
include("metropolis_hastings_trace.jl")
include("rna.jl")
include("hmm.jl")



"""
test(r,transitions,G,nhist,nalleles,onstates,range)

returns dwell time and mRNA histograms for simulations and chemical master equation solutions

G = 2
r = [.01,.02,.03,.04,.001,.001]
transitions = ([1,2],[2,1],[2,3],[3,1])
OFF,ON,mhist,modelOFF,modelON,histF = test(r,transitions,G,20,1,[2],collect(1.:200))

e.g.

julia> incluce("test.jl")

julia> G = 3;

julia> r = [.01,.02,.03,.04,.001,.001];

julia> transitions = ([1,2],[2,1],[2,3],[3,1]);

julia> onstates = [2,3]

julia> OFFsim,ONsim,mRNAsim,OFFchem,ONchem,mRNAchem = test(r,transitions,G,20,1,onstates,collect(1.:200));


"""

"""
r = [0.05,.1,.001,.001];
transitions = ([1,2],[2,1]);
G = 2;
R = 0;
onstates = [2];
trace=make_trace(r,transitions,G,R,onstates,1.,50,10);
interval = 5/3
data = test_data(trace,interval);
model = test_model(r,transitions,G,onstates);
options = test_options(1000);
@time fit,waic = metropolis_hastings(data,model,options);
@time fit,stats,measures = run_mh(data,model,options);

"""


function test(r,transitions,G,nhist,nalleles,onstates,range)
	OFF,ON,mhist = simulator(r,transitions,G,0,0,nhist,nalleles,onstates=onstates,range=range)
	modelOFF,modelON,histF = test_cm(r,transitions,G,nhist,nalleles,onstates,range)
	OFF,ON,mhist,modelOFF,modelON,histF
end

function test_cm(r,transitions,G,nhist,nalleles,onstates,range)
	components=make_components(transitions,G,r,nhist,set_indices(length(transitions)),onstates)
	T = make_mat_T(components.tcomponents,r)
	TA = make_mat_TA(components.tcomponents,r)
	TI = make_mat_TI(components.tcomponents,r)
	M = make_mat_M(components.mcomponents,r)
	modelOFF,modelON = offonPDF(T,TA,TI,range,r,G,transitions,onstates)
	histF = steady_state(M,components.mcomponents.nT,nalleles,nhist)
	modelOFF,modelON,histF
end

test_sim(r,transitions,G,nhist,nalleles,onstates,range) = simulator(r,transitions,G,0,0,nhist,nalleles,onstates=onstates,range=range)

function make_trace(r,transitions,G,R,onstates,interval,steps,ntrials)
	trace = Array{Array{Float64}}(undef,ntrials)
	for i in eachindex(trace)
		trace[i] = simulator(r,transitions,G,R,0,1,1,onstates=onstates,traceinterval=interval,totalsteps=steps)[1:end-1,2]
		for t in eachindex(trace[i])
			if trace[i][t] < .5
				trace[i][t] = 30 + 14*randn()
			else
				trace[i][t] = 230 + 14*randn() + 75*randn()
			end
		end
	end
	trace
end

function test_data(trace,interval)
	TraceData("trace","test",interval,trace)
end

function test_model(r,transitions,G,onstates,propcv=.05,f=Normal,cv=1.)
	components=make_components(transitions,G,r,2,set_indices(length(transitions)),onstates)
	fittedparam=collect(1:length(transitions))
	decayprior = ejectprior = 1.
	nsets = 1
	d = prior_rna(r,G,nsets,fittedparam,decayprior,ejectprior,f,cv)
	GMmodel{typeof(r),typeof(d),typeof(propcv),typeof(fittedparam),typeof(method),typeof(components)}(G,1,r,d,propcv,fittedparam,method,transitions,components,onstates)
end

function test_options(samplesteps::Int=100000,warmupsteps=0,annealsteps=0,maxtime=1000.,temp=1.,tempanneal=100.)

	 MHOptions(samplesteps,warmupsteps,annealsteps,maxtime,temp,tempanneal)

end



# fit,stats,measures = run_mh(data,model,options,nchains);

# function make_components_T(transitions,G,R,r)
# 	elementsT = set_elements_T(transitions,G,R,2)
# 	make_components_TAI(elementsT,G,R,2)
# end

# function mat_trm2(transitions,G,R,r)
# 	components = make_components_T(transitions,G,R,r)
# 	T = make_mat_T(components,r)
# 	TA = make_mat_TA(components,r)
# 	TI = make_mat_TI(components,r)
# 	return T,TA,TI
# end


# function trm2(transitions,G,R,r,bins,nhist)
# 		T,TA,TI = mat_trm2(transitions,G,R,r)
#         modelOFF, modelON = offonPDF(T,TA,TI,bins,r,G-1,R,1)
# 		components = make_components_M(transitions,G,R,nhist,r[2*G+R],2)
#         M = make_mat_M(components,r)
#         histF = steady_state(M,components.nT,1,nhist)
#     	return modelOFF, modelON, histF
# end

# function sims(G,R,r,bins,nhist)
# 	SI,SA = telegraphprefast(bins,G-1,R,r,10000000,10000,1)
# 	hist = telegraphprefast(G-1,R,r,10000000,1e8,nhist,1)
# 	return pdf_from_cdf(bins,SI), pdf_from_cdf(bins,SA),hist
# end


# function mat_GR_DT(r,n,nr)
#     gammap,gamman = get_gamma(r,n)
#     nu = get_nu(r,n,nr)
#     T,B = transition_rate_mat(n,nr,gammap,gamman,nu)
#     T,TA,TI = transition_rate_mat(T,n,nr,2)
# 	return sparse(T),sparse(TA),sparse(TI)
# end
# function trm(G,R,r,bins,nhist)
# 		T,TA,TI = mat_GR_DT(r,G-1,R)
#         modelOFF, modelON = offonPDF(T,TA,TI,bins,r,G-1,R,1)
#         histF = steady_state(r,G-1,R,nhist,1)
#     return modelOFF, modelON, histF
# end
