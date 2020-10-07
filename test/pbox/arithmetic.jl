######
# This file is part of the ProbabilityBoundsAnalysis.jl package.
#
#   Arithmetic tests
#
######


@testset "Arithmetic" begin
    
    mean1 = 2;  mean2 = 2;
    sigma1 =4;  sigma2 = 4;

    a = N(mean1, sigma1); b = N(mean2, sigma2);

    sig(sig1, sig2, cor) = sqrt(sig1^2 + sig2^2 + 2*cor*sig1*sig2)
    meanOut = mean1 + mean2;    

    Nsamples = 100; Ncdf = 100;

    @testset "Independent" begin

        c = a + b;
        d = conv(a,b,op = +,corr = 0);
        e = convIndep(a,b, op = +)

        sigma = sig(sigma1,sigma2,0)
        correct = Normal(meanOut, sigma)

        #@test e.shape == d.shape == "normal"

        @test all(c.u .== d.u .== e.u)
        @test all(c.d .== d.d .== e.d)

        @test mean(d) == mean(e)
        @test mean(correct) ∈ mean(d)
        @test var(d) == var(e)
        @test 32 ∈ var(d)
        
        @test all(map(!,c.bounded))


        Urand = rand(Nsamples,1);
        pbaSamps = cut.(d,Urand);
        DistSamp = quantile.(correct,Urand);

        @test all(DistSamp .∈ pbaSamps)  

    end

    @testset "Perfect" begin

        d = conv(a,b,op = +,corr = 1);
        e = convPerfect(a,b, op = +)

        sigma = sig(sigma1,sigma2,1)
        correct = Normal(meanOut, sigma)

        #@test e.shape == d.shape == "normal"

        @test all(d.u .== e.u)
        @test all(d.d .== e.d)

        @test mean(d) == mean(e)
        @test mean(correct) ∈ mean(d)
        @test var(d) == var(e)
        @test var(correct) ∈ var(d)

        @test all(map(!,d.bounded))


        Urand = rand(Nsamples,1);
        pbaSamps = cut.(d,Urand);
        DistSamp = quantile.(correct,Urand);
        
        @test all(DistSamp .∈ pbaSamps)  

        cdfPoints = range(d.u[1], stop = d.d[end], length = Ncdf);
        pbaCdf = cdf.(d,cdfPoints);
        DistCdf = cdf.(correct,cdfPoints);

        @test all(DistCdf .∈ pbaCdf)  

    end

    @testset "Opposite" begin

        d = conv(a,b, op = +, corr = -1);
        e = convOpposite(a,b, op = +)

        sigma = sig(sigma1,sigma2,-1)
        correct = Normal(meanOut, sigma)

        #@test e.shape == d.shape == "normal"

        @test all(d.u .== e.u)
        @test all(d.d .== e.d)

        @test mean(d) == mean(e)
        @test mean(correct) ∈ mean(d)
        @test var(d) == var(e)
        @test var(correct) ∈ var(d)

        @test all(map(!,d.bounded))


        Urand = rand(Nsamples,1);
        pbaSamps = cut.(d,Urand);
        DistSamp = quantile.(correct,Urand);
        
        @test all(DistSamp .∈ pbaSamps)  

        cdfPoints = range(d.u[1], stop = d.d[end], length = Ncdf);
        pbaCdf = cdf.(d,cdfPoints);
        DistCdf = cdf.(correct,cdfPoints);

        @test all(DistCdf .∈ pbaCdf)  

    end

    @testset "correlated" begin

        R = 0.5
        d = conv(a,b, op = +,corr = R);
        e = sigma(a,b, C = GauCopula(R), op= +)

        sigCor = sig(sigma1,sigma2,R)
        correct = Normal(meanOut, sigCor)

        #@test e.shape == d.shape == "normal"

        @test all(d.u .== e.u)
        @test all(d.d .== e.d)

        @test mean(d) == mean(e)
        @test mean(correct) ∈ mean(d)
        @test var(d) == var(e)
        @test var(correct) ∈ var(d)

        @test all(map(!,d.bounded))


        Urand = rand(Nsamples,1);
        pbaSamps = cut.(d,Urand);
        DistSamp = quantile.(correct,Urand);
        
        @test all(DistSamp .∈ pbaSamps)  

        cdfPoints = range(d.u[1], stop = d.d[end], length = Ncdf);
        pbaCdf = cdf.(d,cdfPoints);
        DistCdf = cdf.(correct,cdfPoints);

        @test all(DistCdf .∈ pbaCdf)  
    end


    @testset "Frechet" begin

        d = conv(a,b, op = +, corr = interval(-1,1));
        e = convFrechet(a,b, op = +)

        #@test e.shape == d.shape == "normal"

        @test all(d.u .== e.u)
        @test all(d.d .== e.d)

        @test mean(d) == mean(e)
        @test var(d) == var(e)

        @test all(map(!,d.bounded))
        
    end

    ###
    #   Need to add testing for:
    #       ->  Shapes
    #       ->  Verify with Monte Carlo?
    #       ->  Double loop Mc for Frechet?
    ###

end