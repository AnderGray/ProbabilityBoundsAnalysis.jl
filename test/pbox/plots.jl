######
# This file is part of the ProbabilityBoundsAnalysis.jl package.
#
#   Plotting tests
#
######


@testset "Plots" begin

    @testset "p-box plots" begin

        a = U(0, 1)

        @test_throws UndefVarError plot(a)

        using PyPlot

        @test @isdefined plot(a, plotting = false)

        @test @isdefined plot(a, name = "same", plotting = false)
        @test @isdefined plot(a, col = "red", plotting = false)
        @test @isdefined plot(a, alpha = 0.5, plotting = false)
        @test @isdefined plot(a, fontsize = 20, plotting = false)
        @test @isdefined plot(a, false, plotting = false)

        b = interval(2, 3)

        @test @isdefined plot(b, plotting = false)
        @test @isdefined plot(b, name = "same", plotting = false)
        @test @isdefined plot(b, col = "red", plotting = false)
        @test @isdefined plot(b, alpha = 0.5, plotting = false)
        @test @isdefined plot(b, fontsize = 20, plotting = false)
        @test @isdefined plot(b, false, plotting = false)

    end

    @testset "copula and bivariate plots" begin

        a = πCop()

        PyPlot.ioff()

        @test @isdefined plot(a)
        @test @isdefined plot(a, name = "same")
        @test @isdefined plot(a, fontsize = 18, alpha = 0.7)
        @test @isdefined plotStep(a)

        b = N(0,1)

        j1 = bivpbox(b, b, a);

        @test @isdefined plot(j1)

        samps = ProbabilityBoundsAnalysis.sample(j1, 10)

        @test @isdefined scatter(samps)
    end

    PyPlot.ion()
end
