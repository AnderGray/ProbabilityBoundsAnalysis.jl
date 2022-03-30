######
# This file is part of the ProbabilityBoundsAnalysis.jl package.
#
#   Distribution tests
#
######

using Distributions

@testset "Distributions" begin

    Nsamples = 1000; Ncdf = 1000;

    function test_dist(x :: pbox, a)

        samps = rand(Nsamples);

        @test length(x.d) == ProbabilityBoundsAnalysis.parametersPBA.steps

        pbaSamples = cut.(x,samps);
        distSamps  = quantile.(a,samps);

        @test all( distSamps .∈ pbaSamples )

        testVals = range(left(x), stop = right(x),length = Ncdf)
        cdfPba = cdf.(x, testVals);
        cdfDis = cdf.(a, testVals);

        @test all( cdfDis .∈ cdfPba)

        EndPoints = rand(Nsamples,2) .* (right(x) .- left(x)) .+ left(x);
        EndPoints = sort(EndPoints, dims = 2);

        massInts = interval.(EndPoints[:,1], EndPoints[:,2])
        massPba = mass.(x,massInts)
        massDis = cdf.(a,EndPoints[:,2]) .- cdf.(a, EndPoints[:,1])

        @test all( massDis .∈ massPba)

    end

    @testset "Set based operations" begin

        a = U(interval(0, 1), 2);
        b = U(1, 2);

        c = imp(a, b)

        @test c.u == b.u
        @test c.d == b.d
        @test mean(c) == mean(b)
        @test var(c) == var(b)

        b2 = U(1, interval(2, 3))

        c2 = imp(a, b2)

        @test c2.u == b.u
        @test c2.d == b.d
        @test mean(c2) == mean(b)
        @test ~isempty(var(c2) ∩ var(b))

        b3 = U(2, 3)

        @test_throws ArgumentError imp(a, b3)

        b = U(2, 3)
        c = env(a, b)

        @test c.d == b.d
        @test c.u == a.u
        @test mean(c) == hull(mean(a), mean(b))
        @test var(c) == hull(var(a), var(b))

        b = makepbox(1..2)
        c = env(a, b)

        @test c.d == b.d
        @test c.u[1:100] == a.u[1:100]
        @test c.u[101:end] == b.u[101:end]
        @test mean(c) == hull(mean(a), mean(b))
        @test var(c) == hull(var(a), var(b))

    end

    @testset "Parametric p-boxes - 200 steps" begin

        test_dist(normal(0,1), Normal(0,1))
        test_dist(normal(0..1, 1..2), Normal(0,1))

        test_dist(uniform(0, 1), Uniform(0,1))
        test_dist(uniform(0..0.5, 1..1), Uniform(0,1))

        test_dist(beta(2,3), Beta(2, 3))
        test_dist(beta(2..3, 3..4), Beta(2, 3))

        test_dist(betaPrime(2, 3), BetaPrime(2, 3))
        test_dist(betaPrime(2..3, 3..4), BetaPrime(2, 3))

        test_dist(biweight(1,2), Biweight(1, 2))
        test_dist(biweight(1..2,2..2), Biweight(1, 2))

        test_dist(cauchy(0, 1), Cauchy(0, 1))
        test_dist(cauchy(0..0.5, 1), Cauchy(0, 1))

        test_dist(chi(1), Chi(1))
        test_dist(chi(1..2), Chi(1))

        test_dist(chisq(1), Chisq(1))
        test_dist(chisq(1..2), Chisq(1))

        test_dist(cosine(1,2), Cosine(1,2))
        test_dist(cosine(1..2,2), Cosine(1,2))

        test_dist(epanechnikov(1,2), Epanechnikov(1,2))
        test_dist(epanechnikov(1..2,2), Epanechnikov(1,2))

        test_dist(exponential(1), Exponential(1))
        test_dist(exponential(1..2), Exponential(1))

        test_dist(fDist(2, 3), FDist(2, 3))
        test_dist(fDist(2..3, 3), FDist(2, 3))

        test_dist(gamma(1, 1), Gamma(1, 1))
        test_dist(gamma(1, 1..2), Gamma(1, 1))

        test_dist(laplace(0, 1), Laplace(0, 1))
        test_dist(laplace(0..1, 1), Laplace(0, 1))

        test_dist(levy(0, 1), Levy(0, 1))
        test_dist(levy(0..1, 1), Levy(0, 1))

        test_dist(lognormal(2,1), ProbabilityBoundsAnalysis.pbaLogNormal(2,1))
        test_dist(lognormal(2..3,1), ProbabilityBoundsAnalysis.pbaLogNormal(2,1))

        ### Broken
        #test_dist(frechet(1, 1), Frechet(1, 1))
        #test_dist(frechet(1..2, 1), Frechet(1, 1))
        #test_dist(ksdist(2), KSDist(2))

        #test_dist(erlang(1,1), Erlang(1,1))
        #test_dist(erlang(1..2,1), Erlang(1,1))
    end

    @testset "Parametric p-boxes - 20 steps" begin

        ProbabilityBoundsAnalysis.setSteps(20)

        test_dist(normal(0,1), Normal(0,1))
        test_dist(normal(0..1, 1..2), Normal(0,1))

        test_dist(uniform(0, 1), Uniform(0,1))
        test_dist(uniform(0..0.5, 1..1), Uniform(0,1))

        test_dist(beta(2,3), Beta(2, 3))
        test_dist(beta(2..3, 3..4), Beta(2, 3))

        test_dist(betaPrime(2, 3), BetaPrime(2, 3))
        test_dist(betaPrime(2..3, 3..4), BetaPrime(2, 3))

        test_dist(biweight(1,2), Biweight(1, 2))
        test_dist(biweight(1..2,2..2), Biweight(1, 2))

        test_dist(cauchy(0, 1), Cauchy(0, 1))
        test_dist(cauchy(0..0.5, 1), Cauchy(0, 1))

        test_dist(chi(1), Chi(1))
        test_dist(chi(1..2), Chi(1))

        test_dist(chisq(1), Chisq(1))
        test_dist(chisq(1..2), Chisq(1))

        test_dist(cosine(1,2), Cosine(1,2))
        test_dist(cosine(1..2,2), Cosine(1,2))

        test_dist(epanechnikov(1,2), Epanechnikov(1,2))
        test_dist(epanechnikov(1..2,2), Epanechnikov(1,2))

        test_dist(exponential(1), Exponential(1))
        test_dist(exponential(1..2), Exponential(1))

        test_dist(fDist(2, 3), FDist(2, 3))
        test_dist(fDist(2..3, 3), FDist(2, 3))

        test_dist(gamma(1, 1), Gamma(1, 1))
        test_dist(gamma(1, 1..2), Gamma(1, 1))

        test_dist(laplace(0, 1), Laplace(0, 1))
        test_dist(laplace(0..1, 1), Laplace(0, 1))

        test_dist(levy(0, 1), Levy(0, 1))
        test_dist(levy(0..1, 1), Levy(0, 1))

        test_dist(lognormal(2,1), ProbabilityBoundsAnalysis.pbaLogNormal(2,1))
        test_dist(lognormal(2..3,1), ProbabilityBoundsAnalysis.pbaLogNormal(2,1))

        ProbabilityBoundsAnalysis.setSteps(200)

    end


    @testset "c-boxes" begin

        @test_throws ArgumentError KN(-1, 2)
        @test_throws ArgumentError KN(3, 2)

        @test all(KN(0, 0).d .== 1)
        @test all(KN(0, 0).u .== 0)

        @test all(KN(0, 1).u .== 0)
        @test all(KN(1, 1).d .== 1)

        @test all(KN(0, 1..2).u .== 0)
        @test all(KN(0..1, 1).d .== 1)

    end

    @testset "Distribution-free p-boxes" begin



    end

end
