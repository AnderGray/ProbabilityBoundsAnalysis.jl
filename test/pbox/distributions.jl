######
# This file is part of the ProbabilityBoundsAnalysis.jl package.
#
#   Distribution tests
#
######

using Distributions

@testset "Distributions" begin

    Nsamples = 5000; Ncdf = 5000;

    @testset "Normal" begin

        x = normal(interval(0,1),interval(1,2));    

        @test mean(x) == interval(0,1);
        @test var(x) == interval(1,2)^2;
        @test x.shape == "normal"
        @test all( map(!,x.bounded))

        x = normal(0,1);
        a = Normal(0,1);

        samps = rand(Nsamples);

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

    @testset "Uniform" begin

    x = U(interval(0,1),interval(1,2));    

    @test x.u[1] == 0;
    @test x.d[end] == 2;
    @test x.shape == "uniform"
    @test all(x.bounded)

    x = uniform(0,1);
    a = Uniform(0,1);

    samps = rand(Nsamples);

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

    @testset "Beta" begin

    x = beta(interval(1,3),interval(1,2));    

    @test x.shape == "beta"
    
    x = uniform(0,1);
    a = Uniform(0,1);

    samps = rand(Nsamples);

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

    # Distribution Free constructors need alot of test
    #


end