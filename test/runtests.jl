using StochasticGene
using Test

function fit_rna_test()

    data = make_data("testgene","WT","data/",false,"test",".")
    model = make_model("testgene","testcell",2,false,[1,2,3],(),"test","test",1,".",data)
    options = StochasticGene.MHOptions(1000,100,1000,1000,60.,1.,100.)
    param = StochasticGene.get_param(model)
    fit,stats,waic = StochasticGene.run_mh(data,model,options,nchains);
    nothing
end


function teststeadystatemodel(model::AbstractGMmodel,nhist)
    G = model.G
    r = model.rates
    g1 = steady_state(r[1:2*G],G-1,nhist,model.nalleles)
    g2 = simulatorGM(r[1:2*G],G-1,nhist,model.nalleles)
    
    return g1,g2
end

@testset "StochasticGene" begin
    # Write your tests here.
end
