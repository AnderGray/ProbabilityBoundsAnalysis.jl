######
# This file is part of the ProbabilityBoundsAnalysis.jl package.
#
#   Copula tests
#
######

using ProbabilityBoundsAnalysis: Frechet, massD, massU

Ngrid = 1000;
πfunc(x,y) = x*y
mfunc(x,y) = min(x, y)
wfunc(x,y) = max(x + y - 1, 0)

function test_copula(C, func)

    xs = range(0, 1, length = Ngrid)

    @test func(0.4, 0.4) ∈ cdf(C,0.4,0.4)
    @test func(0.4, 0.4) ∈ C(0.4,0.4)

    cdfs = cdf(C, xs, xs)
    cdfs_lo = cdfD(C, xs, xs)
    cdfs_hi = cdfU(C, xs, xs)

    real_vals = [func(x, y) for x in xs, y in xs]

    @test all(real_vals .∈ cdfs)
    @test all(real_vals .>= getfield.(cdfs_lo, :lo))
    @test all(real_vals .<=  getfield.(cdfs_hi, :hi))

    @test all(func.(0.4, xs) .∈ C(0.4,xs))
    @test all(func.(xs, 0.4) .∈ C(xs, 0.4))

    @test all(func.(0.4, xs) .∈ C(0.4,collect(xs)))
    @test all(func.(collect(xs), 0.4) .∈ C(xs, 0.4))

    @test all(func.(xs, xs) .∈ C.(xs, xs))

    x_ints = mince(0..1, 100)

    @test all(range(0,1, length = 100) .∈ C(1, x_ints))
    @test all(range(0,1, length = 100) .∈ C(x_ints, 1))

    @test all(zeros(100) .∈ C(0, x_ints))
    @test all(zeros(100) .∈ C(x_ints, 0))

    cdf_ints = [C(x, y) for x in x_ints, y in x_ints]
    cdf_int_real = [func(x, y) for x in x_ints, y in x_ints]

    @test all(cdf_int_real .⊆ cdf_ints)

    cop_masses = [mass(C, interval(0, x), interval(0, y)) for x in xs, y in xs]

    cop_massesU = [massU(C, interval(0, x), interval(0, y)) for x in xs, y in xs]
    cop_massesD = [massD(C, interval(0, x), interval(0, y)) for x in xs, y in xs]

    @test all(real_vals .∈ cop_masses)
    @test all(cop_massesU .⊆ cop_masses)
    @test all(cop_massesD .⊆ cop_masses)
end

@testset "Copulas" begin

    @testset "π" begin
        a = πCop()
        test_copula(a, πfunc)
    end

    @testset "W" begin
        a = W()
        test_copula(a, wfunc)
    end

    @testset "M" begin
        a = M()
        test_copula(a, mfunc)
    end

    @testset "Frank" begin
        a = Frank(0)
        test_copula(a, πfunc)

        a = Frank(-Inf)
        test_copula(a, wfunc)

        a = Frank(Inf)
        test_copula(a, mfunc)
    end

    @testset "Frank - interval" begin
        a = Frank(-1 .. 1)
        test_copula(a, πfunc)
    end

    @testset "Clayton" begin
        a = Clayton(0)
        test_copula(a, πfunc)

        a = Clayton(-1)
        test_copula(a, wfunc)

        a = Clayton(Inf)
        test_copula(a, mfunc)
    end

    @testset "Clayton - interval" begin
        a = Clayton(0 .. 1)
        test_copula(a, πfunc)
    end

    @testset "Gaussian " begin

        a = GauCopula(0)
        test_copula(a, πfunc)

        a = GauCopula(1)
        test_copula(a, mfunc)

        a = GauCopula(-1)
        test_copula(a, wfunc)

    end

    @testset "Gaussian - interval" begin

        a = GauCopula(interval(0,1))
        test_copula(a, πfunc)
        test_copula(a, mfunc)

        a = GauCopula(interval(-1,0))
        test_copula(a, πfunc)
        test_copula(a, wfunc)

    end

    @testset "Frechet bounds" begin

        a = env(W(), M())
        test_copula(a, πfunc)
        test_copula(a, mfunc)
        test_copula(a, wfunc)

        b = Frechet()
        test_copula(b, πfunc)
        test_copula(b, mfunc)
        test_copula(b, wfunc)

    end

    @testset "conditionals" begin

        U1 = U(0,1);
        a = πCop()
        xs = range(0.1, 0.9, length=10)

        a_condX = conditionalX.(a, xs);
        a_condY = conditionalY.(a, xs);

        @test all(imp.(U1, a_condX) .⊆ U1)
        @test all(imp.(U1, a_condY) .⊆ U1)

        xs = range(0.1, 0.9, length=10)
        m = M()

        m_condX = conditionalX.(m, xs);
        m_condY = conditionalY.(m, xs);

        ints = xs .± 0.1

        @test all( 1 .∈ mass.(m_condX, ints))
        @test all( 1 .∈ mass.(m_condY, ints))

        w = W()
        w_condX = conditionalX.(w, xs);
        w_condY = conditionalY.(w, xs);

        ints = (1 .- xs) .± 0.1

        @test all( 1 .∈ mass.(w_condX, ints))
        @test all( 1 .∈ mass.(w_condY, ints))

    end

end
