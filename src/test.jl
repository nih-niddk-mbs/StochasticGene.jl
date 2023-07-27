
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



"""
test_cm(r,transitions,G,nhist,nalleles,onstates,range)

G = 2
r = [.01,.02,.03,.04,.001,.001]
transitions = ([1,2],[2,1],[2,3],[3,1])
OFF,ON,mhist,modelOFF,modelON,histF = test(r,transitions,G,20,1,[2],collect(1.:200))

e.g.

julia> G = 3;

julia> r = [.01,.02,.03,.04,.001,.001];

julia> transitions = ([1,2],[2,1],[2,3],[3,1]);

julia> OFF,ON,mhist,modelOFF,modelON,histF = test(r,transitions,G,20,1,[2],collect(1.:200))
Vector{Float64}
([0.004986113758270482, 0.005144376072099545, 0.005229463337599042, 0.005384322160808125, 0.005770618346175838, 0.005894845753805103, 0.006194352928363329, 0.0062028616549132795, 0.006328790807852534, 0.006464930432651728  …  0.002147602581207286, 0.002169725270237155, 0.0021424973452773164, 0.0021322868734173768, 0.001968919323658344, 0.002103357203147548, 0.0019110599831186866, 0.0020267786641980016, 0.0019927437579982028, 0.0019570071064884146], [0.04907510223979633, 0.0465919830315436, 0.0441145655895555, 0.04160293754997955, 0.04008484228201791, 0.0378982149195267, 0.03618911048169947, 0.034522769290857076, 0.03263263377412738, 0.03124140280555409  …  5.70176626464463e-6, 8.552649396966945e-6, 1.4254415661611575e-6, 1.4254415661611575e-6, 0.0, 1.4254415661611575e-6, 2.850883132322315e-6, 0.0, 0.0, 1.4254415661611575e-6], [0.8930025545901369, 0.10070712078313866, 0.0060716276315020395, 0.00021869699522261148, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.011507484786359158, 0.011392901825339419, 0.011279679874098754, 0.011167140664500235, 0.011056414072153693, 0.010945649502079696, 0.010836171264803546, 0.010727857953268003, 0.010620709567473097, 0.010514726107418765  …  0.0017214329316598287, 0.0017043436375838304, 0.0016874167250036287, 0.0016706497898014575, 0.0016540404278600647, 0.001637586235061813, 0.0016212848072893217, 0.001605133740424954, 0.001589130630351457, 0.001573273072951194], [0.020108643762429337, 0.02025617771086947, 0.02034867780719349, 0.020395098487701244, 0.020398017630536794, 0.02036288587463498, 0.02029093028656673, 0.02018679014884625, 0.020052338748856313, 0.01989095489498634  …  5.155985894233495e-5, 4.961612868379767e-5, 4.774434639820735e-5, 4.594088667948768e-5, 4.420212412078449e-5, 4.2524433315799205e-5, 4.090418885778878e-5, 3.933811015378926e-5, 3.78419504803493e-5, 3.641097251853732e-5], [0.8957142251674732, 0.09778204909443922, 0.0061948505211017624, 0.0002965844677297766, 1.1860253802978403e-5, 4.168286955594493e-7, 1.3265972845065816e-8, 3.8978275758618964e-10, 1.0714734926772281e-11, 2.7845689471113544e-13, 7.11684540581175e-15, 3.8566612683656346e-16, 2.0819810144896226e-16, 1.887009673667293e-16, 1.7502292428581908e-16, 1.6325359423904853e-16, 1.5297439726137662e-16, 1.4393099534679304e-16, 1.3590886508584802e-16, 1.2650622706296916e-16])


"""



function test_cm(r,transitions,G,nhist,nalleles,onstates,range)
	components=make_components(transitions,G,r,nhist,set_indices(length(transitions)),onstates)
	T = make_mat_T(components.tcomponents,r)
	TA = make_mat_TA(components.tcomponents,r)
	TI = make_mat_TI(components.tcomponents,r)
	M = make_mat_M(components.mcomponents,r)
	offonPDF(TA,TI,range,G,transitions,onstates)
	histF = steady_state(M,components.mcomponents.nT,nalleles,nhist)
	modelOFF,modelON,histF
end

test_sim(r,transitions,G,nhist,nalleles,onstates,range) = simulator(r,transitions,G,0,0,nhist,nalleles,onstates=onstates,range=range)






# function test(r,transitions,G,nhist,nalleles,onstates,range)
# 	OFF,ON,mhist = simulator(r,transitions,G,0,0,nhist,nalleles,onstates=onstates,range=range)
# 	modelOFF,modelON,histF = test_cm(r,transitions,G,nhist,nalleles,onstates,range)
# 	OFF,ON,mhist,modelOFF,modelON,histF
# end




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
