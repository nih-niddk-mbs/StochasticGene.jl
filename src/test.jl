include("/Users/carsonc/github/StochasticGene/src/common.jl")
include("/Users/carsonc/github/StochasticGene/src/transition_rate_matrices2.jl")
include("/Users/carsonc/github/StochasticGene/src/transition_rate_matrices.jl")
include("/Users/carsonc/github/StochasticGene/src/chemical_master2.jl")
include("/Users/carsonc/github/StochasticGene/src/GillespieSimulator.jl")
include("/Users/carsonc/github/StochasticGene/src/utilities.jl")

using DifferentialEquations
using LSODA


function make_components_T(transitions,G,R,r)
	elementsT = set_elements_T(transitions,G,R,2)
	make_components_TAI(elementsT,G,R,2)
end

function mat_trm2(transitions,G,R,r)
	components = make_components_T(transitions,G,R,r)
	T = make_mat_T(components,r)
	TA = make_mat_TA(components,r)
	TI = make_mat_TI(components,r)
	return T,TA,TI
end


function trm2(transitions,G,R,r,bins,nhist)
		T,TA,TI = mat_trm2(transitions,G,R,r)
        modelOFF, modelON = offonPDF(T,TA,TI,bins,r,G-1,R,1)
		components = make_components_M(transitions,G,R,nhist,r[2*G+R],2)
        M = make_mat_M(components,r)
        histF = steady_state(M,components.nT,1,nhist)
    	return modelOFF, modelON, histF
end

function sims(G,R,r,bins,nhist)
	SI,SA = telegraphprefast(bins,G-1,R,r,10000000,10000,1)
	hist = telegraphprefast(G-1,R,r,10000000,1e8,nhist,1)
	return pdf_from_cdf(bins,SI), pdf_from_cdf(bins,SA),hist
end


function mat_GR_DT(r,n,nr)
    gammap,gamman = get_gamma(r,n)
    nu = get_nu(r,n,nr)
    T,B = transition_rate_mat(n,nr,gammap,gamman,nu)
    T,TA,TI = transition_rate_mat(T,n,nr,2)
	return sparse(T),sparse(TA),sparse(TI)
end
function trm(G,R,r,bins,nhist)
		T,TA,TI = mat_GR_DT(r,G-1,R)
        modelOFF, modelON = offonPDF(T,TA,TI,bins,r,G-1,R,1)
        histF = steady_state(r,G-1,R,nhist,1)
    return modelOFF, modelON, histF
end
