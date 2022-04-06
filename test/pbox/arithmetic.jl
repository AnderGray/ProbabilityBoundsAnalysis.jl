######
# This file is part of the ProbabilityBoundsAnalysis.jl package.
#
#   Arithmetic tests
#
######


@testset "Arithmetic" begin

    mean1 = 2;  mean2 = 2;
    sigma1 = 4;  sigma2 = 4;

    a = N(mean1, sigma1); b = N(mean2, sigma2);

    sig(sig1, sig2, cor) = sqrt(sig1^2 + sig2^2 + 2*cor*sig1*sig2)
    meanOut = mean1 + mean2;

    Nsamples = 5000; Ncdf = 5000;

    @testset "Independent" begin

        c = conv(a,b,op = +,corr = 0);
        d = convIndep(a,b, op = +)

        sigma = sig(sigma1, sigma2, 0)
        correct = Normal(meanOut, sigma)

        #@test e.shape == d.shape == "normal"

        @test all(c.u .== d.u)
        @test all(c.d .== d.d)

        @test mean(c) == mean(d)
        @test var(c) == var(d)
        @test mean(correct) ∈ mean(d)
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

        ops = [+, -, *, /, min, max];

        for op in ops

            x1 = U(-1,1)
            x2 = U(1,2)

            c = convFrechet(x1, x2, op = op)

            cI = convIndep(x1, x2, op = op)
            cP = convPerfect(x1, x2, op = op)
            cO = convOpposite(x1, x2, op = op)

            @test all( c.u .<= cI.u) && all( c.d .>= cI.d)
            @test all( c.u .<= cP.u) && all( c.d .>= cP.d)
            @test all( c.u .<= cO.u) && all( c.d .>= cO.d)
        end
    end

    @testset "Partially known dependence" begin

        x1 = U(1, 3)
        x2 = N(10, 1)

        ops = [+, *, min, max];

        Cs = Vector{copula}(undef, 9)

        Cs[1] = M()
        Cs[5] = πCop()
        Cs[end] = W()

        corrs = [0.8, 0.5, 0.2, -0.2, -0.5, -0.8]
        Gs = GauCopula.(corrs)

        Cs[2:4] = Gs[1:3]
        Cs[6:8] = Gs[4:end]

        for op in ops

            partial = [tauRho(x1, x2, op = op, C = C) for C in Cs]
            precise = [sigma(x1, x2, op = op, C = C) for C in Cs]

            precise[1] = convPerfect(x1, x2, op = op)
            precise[end] = convOpposite(x1, x2, op = op)

            for i = 1:length(partial)
                for j = 1:i
                    @test imp(precise[j], partial[i]) ⊆ partial[i]
                    @test partial[j] ⊆ partial[i]
                end
            end
        end
    end

    @testset "P-box unary" begin

        a = U(1,2)
        aneg = -a
        aRecep = 1/a

        c0 = convOpposite(a, aneg, op=+)
        c1 = convOpposite(a, aRecep, op=*)

        @test  c0(0) == interval(0,1)
        @test  c1(1) == interval(0,1)

    end
end
