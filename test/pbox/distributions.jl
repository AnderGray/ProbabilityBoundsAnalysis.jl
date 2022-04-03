######
# This file is part of the ProbabilityBoundsAnalysis.jl package.
#
#   Distribution tests
#
######

using Distributions

@testset "Distributions" begin

    Nsamples = 101; Ncdf = 101;

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

        x1 = uniform(mean = 0, std = 1..2);
        @test mean(x1) == 0;
        @test 1..2 ⊆ std(x1);

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

    @testset "Moment constraints" begin

        moms = checkMomentsAndRanges(0, 2)

        @test moms[1].lo == 0  # min
        @test moms[2].hi == 2  # max
        @test moms[3] == 0..2  # mean
        @test moms[4] == 0..1  # var

        moms = checkMomentsAndRanges(0, 2, 0.5)

        @test moms[3] == 0.5  # mean
        @test moms[4] == 0..0.75  # var

        moms = checkMomentsAndRanges(0, 2, 0.5, 0.5)

        @test moms[3] == 0.5  # mean
        @test moms[4] == 0.5  # var

        @test_throws ArgumentError checkMomentsAndRanges(2, 0)
        @test_throws ArgumentError checkMomentsAndRanges(0, 2, interval(3, 4))
        @test_throws ArgumentError checkMomentsAndRanges(0, 2, 1, 10)


        # checks if mean is reduced
        @test checkMomentsAndRanges(0, 2, interval(1, 4))[3] == interval(1, 2)

        # checks if var is reduced
        @test checkMomentsAndRanges(0, 2, 0.5, interval(0.5, 10))[4] == interval(0.5, 0.75)
    end

    @testset "Distribution-free p-boxes" begin

        x1 = meanVar(0, 1)
        @test mean(x1) == 0
        @test var(x1) == 1

        x2 = meanVar(-0.5..0.5, 1..2)
        @test mean(x2) == interval(-0.5, 0.5)
        @test var(x2) == interval(1..2)

        @test x1 ⊆ x2

        x3 = meanMin(3, 0)

        @test mean(x3) == 3
        @test x3.u[1] == 0

        x4 = meanMin(interval(2,4), 0..1)

        @test mean(x4) == 2..4
        @test x4.u[1] == 0

        @test x3 ⊆ x4

        x5 = meanMinMax(3, 0, 5)
        @test x5 ⊆ x3

        x6 = meanMinMax(2..4, 0, 5)
        @test x5 ⊆ x6
        @test x6 ⊆ x4

        x7 = minMeanVar(0, 3, 0.5)
        @test x7 ⊆ x3

        x8 = minMeanVar(0, 2..4, 0.5..1)
        @test x7 ⊆ x8
        @test x8 ⊆ x4

        x9 = minMaxMeanVar(0, 6, 3, 0.5)
        x10 = meanMax(3, 6)
        @test x9 ⊆ x8
        @test x9 ⊆ x4
        @test x9 ⊆ x10

        x11 = maxMeanVar(10, 0, 0.9..1)
        x12 = uniform(mean = 0, std = 1..1)
        @test x12 ⊆ x11

        @test Chebyshev(0,1) == meanVar(0,1)
        @test chebyshev(0,1) == meanVar(0,1)
        @test cheb(0,1) == MeanVar(0,1)
        @test cheb(0,1) == meanvar(0,1)
        @test cheb(0,4) == meanStd(0,2)

        @test Markov(1, 0) == meanMin(1, 0)
        @test markov(1, 0) == MeanMin(1, 0)
        @test markov(1, 0) == meanmin(1, 0)

        @test meanMax(0, 1) == MeanMax(0, 1)
        @test meanmax(0, 1) == MeanMax(0, 1)

        @test Cantelli(0.5, 0, 1) == meanMinMax(0.5, 0, 1)
        @test cantelli(0.5, 0, 1) == meanminmax(0.5, 0, 1)

        @test mmmv(0, 10, 5, 4) == mmms(0, 10, 5, 2)
        @test minMaxMeanVar(0, 10, 5, 4) == mmmv(0, 10, 5, 4)
        @test Ferson(0, 10, 5, 4) == mmmv(0, 10, 5, 4)
        @test MinMaxMeanVar(0, 10, 5, 4) == mmmv(0, 10, 5, 4)
        @test ferson(0, 10, 5, 4) == mmmv(0, 10, 5, 4)
        @test minmaxmeanvar(0, 10, 5, 4) == minMaxMeanVar(0, 10, 5, 4)
        @test minMaxMeanStd(0, 10, 5, 2) == minMaxMeanVar(0, 10, 5, 4)

    end
end
