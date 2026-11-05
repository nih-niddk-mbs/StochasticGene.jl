using Test
using StochasticGene

@testset "quick-start API" begin
    @testset "rna_histogram" begin
        @test rna_histogram([0, 0, 1, 3]) == [2, 1, 0, 1]
        @test_throws ArgumentError rna_histogram(Int[])
        @test_throws ArgumentError rna_histogram([-1, 2])
    end

    @testset "input validation" begin
        dir = mktempdir()
        @test_throws ArgumentError fit_rna(joinpath(dir, "missing.txt"); outdir=dir)
        @test_throws ArgumentError fit_rna([1, 2, 3, 4]; outdir=dir, halflife=-1.0)
        @test_throws ArgumentError fit_traces([rand(20)]; R=0, outdir=dir)
        @test_throws ArgumentError fit_traces([rand(20)]; interval=0.0, outdir=dir)
        @test_throws ArgumentError fit_traces([rand(3)]; outdir=dir)
        @test_throws ArgumentError fit_traces(joinpath(dir, "no_such_folder"); outdir=dir)
    end

    @testset "trace staging" begin
        dir = mktempdir()
        tracedir = joinpath(dir, "traces")
        n = StochasticGene._stage_traces([rand(20), rand(25)], tracedir, "g", 2.0)
        @test n == 2
        files = sort(readdir(tracedir))
        @test files == ["0001_g_DATA.txt", "0002_g_DATA.txt"]
        t = StochasticGene.readfile(joinpath(tracedir, files[1]))
        @test size(t) == (20, 3)
        @test t[:, 2] ≈ (0:19) .* 2.0
        # restaging replaces earlier files instead of accumulating
        @test StochasticGene._stage_traces([rand(20)], tracedir, "g", 1.0) == 1
        @test length(readdir(tracedir)) == 1
    end

    @testset "fit_rna smoke" begin
        dir = mktempdir()
        h = rna_histogram(rand(0:20, 300))
        r = fit_rna(h; G=2, name="x", halflife=300.0, outdir=dir, samplesteps=2000, nchains=1)
        @test r.rates.parameter == ["Rate12", "Rate21", "Eject", "Decay"]
        @test all(isfinite, r.rates.median)
        @test fit_summary(r.resultdir).median == r.rates.median
    end
end
