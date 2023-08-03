using StochasticGene
using Test

function fit_rna_test(root=".")
    G = 2
    gene = "CENPL"
    cell = "HCT116"
    fish = false
    nalleles = 2
    nsets = 1
    propcv = 0.05
    fittedparam = [1,2,3]
    fixedeffects = ()
    transitions = ([1,2],[2,1])
    ejectprior = 0.05
    r = [0.01, 0.1, 1.0, 0.01006327034802035]
    decayrate = 0.01006327034802035
    datacond = "MOCK"
    datafolder = "data/HCT116_testdata"
    label = "scRNA_test"
    data = StochasticGene.data_rna(gene,datacond,datafolder,fish,label)
    model = StochasticGene.model_rna(data,r,G,nalleles,nsets,propcv,fittedparam,fixedeffects,transitions,decayrate,ejectprior)
    options = StochasticGene.MHOptions(100000,0,0,30.,1.,100.)
    fit,stats,measures = StochasticGene.run_mh(data,model,options,1);
    return stats.meanparam, fit.llml, model
end

function teststeadystatemodel(r,transitions,G,nhist,nalleles)
    g1 = StochasticGene.steady_state(r,transitions,G,nhist,nalleles)
    g2 = StochasticGene.simulator(r,transitions,G,0,0,nhist,nalleles)
    return g1, g2
end


@testset "StochasticGene" begin

        p,ll,model =  fit_rna_test()
        @test isapprox(ll,1766,rtol=0.05)

        r = [0.0014, 0.005, 0.0016, 0.01]
        transitions = ([1,2],[2,1])
        h1,h2 = teststeadystatemodel(r,transitions,2,20,2)

        @test isapprox(h1,h2,rtol=0.05)

end
