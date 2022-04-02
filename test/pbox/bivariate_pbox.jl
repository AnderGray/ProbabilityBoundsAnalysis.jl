######
# This file is part of the ProbabilityBoundsAnalysis.jl package.
#
#   Bivaraite p-box tests
#
######

using ProbabilityBoundsAnalysis: Frechet, massD, massU

Ngrid = 101;
πfunc(x,y) = x*y
mfunc(x,y) = min(x, y)
wfunc(x,y) = max(x + y - 1, 0)

function test_2pbox(J, func)

    xRange = interval(left(J.marg1), right(J.marg1))
    yRange = interval(left(J.marg2), right(J.marg2))

    xs = range(xRange.lo, xRange.hi, length = Ngrid)
    ys = range(yRange.lo, yRange.hi, length = Ngrid)

    @test func(0.4, 0.4) ∈ cdf(J,0.4,0.4)
    @test func(0.4, 0.4) ∈ J(0.4,0.4)

    cdfs = cdf(J, xs, ys)

    real_vals = abs.([func(x, y) for x in xs, y in xs])

    @test all(real_vals .∈ cdfs)

    @test all(func.(0.4, xs) .∈ J(0.4,xs))
    @test all(func.(xs, 0.4) .∈ J(xs, 0.4))

    @test all(func.(0.4, xs) .∈ J(0.4,collect(xs)))
    @test all(func.(collect(xs), 0.4) .∈ J(xs, 0.4))

    @test all(abs.(func.(xs, xs)) .∈ J.(xs, xs))

    x_ints = mince(xRange, Ngrid)
    y_ints = mince(yRange, Ngrid)

    @test all(cdf.(J.marg1, xs) .⊆ J(Inf, x_ints))
    @test all(cdf.(J.marg2, ys) .⊆ J(y_ints, Inf))

    @test all(zeros(Ngrid) .∈ J(-Inf, x_ints))
    @test all(zeros(Ngrid) .∈ J(y_ints, -Inf))

    cdf_ints = [J(x, y) for x in x_ints, y in y_ints]

    @test all(real_vals .⊆ cdf_ints)

    cop_masses = [mass(J, interval(-Inf, x), interval(-Inf, y)) for x in xs, y in xs]

    @test all(real_vals .∈ cop_masses)
end

@testset "bivariate p-boxes" begin

    @testset "Uniform marginals" begin

        x1 = U(0,1); x2 = U(0,1);
        c = πCop()
        j1 = c(x1, x2)

        test_2pbox(j1, πfunc)

        c = M();
        j2 = c(x1,x2)
        test_2pbox(j2, mfunc)

        j3 = W()(x1,x2)
        test_2pbox(j3, wfunc)

    end

    @testset "Normal marginals" begin

        x1 = N(0, 1); x2 = N(0, 1);
        c = πCop();
        j1 = c(x1,x2)

        func(x,y) = ProbabilityBoundsAnalysis.bivariate_cdf(x, y, 0)
        test_2pbox(j1, func)

        c = GauCopula(0.5);
        j2 = c(x1,x2);

        func2(x,y) = ProbabilityBoundsAnalysis.bivariate_cdf(x, y, 0.5)
        test_2pbox(j2, func2)

        c = GauCopula(-0.5);
        j3 = c(x1,x2);

        func3(x,y) = ProbabilityBoundsAnalysis.bivariate_cdf(x, y, -0.5)
        test_2pbox(j3, func3)

        c = Frechet();
        j4 = c(x1,x2);

        func4(x,y) = ProbabilityBoundsAnalysis.bivariate_cdf(x, y, -0.8)
        test_2pbox(j4, func4)

        func5(x,y) = ProbabilityBoundsAnalysis.bivariate_cdf(x, y, -0.3)
        test_2pbox(j4, func5)

        func6(x,y) = ProbabilityBoundsAnalysis.bivariate_cdf(x, y, 0.3)
        test_2pbox(j4, func6)

        func7(x,y) = ProbabilityBoundsAnalysis.bivariate_cdf(x, y, 0.8)
        test_2pbox(j4, func7)

    end

end
