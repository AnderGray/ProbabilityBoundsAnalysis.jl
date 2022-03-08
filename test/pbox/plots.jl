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

        @test plot(a, plotting = false) == nothing

        @test plot(a, name = "same", plotting = false) == nothing
        @test plot(a, col = "red", plotting = false) == nothing
        @test plot(a, alpha = 0.5, plotting = false) == nothing
        @test plot(a, fontsize = 20, plotting = false) == nothing
        @test plot(a, false, plotting = false) == nothing

        b = interval(2, 3)

        @test plot(b, plotting = false) == nothing
        @test plot(b, name = "same", plotting = false) == nothing
        @test plot(b, col = "red", plotting = false) == nothing
        @test plot(b, alpha = 0.5, plotting = false) == nothing
        @test plot(b, fontsize = 20, plotting = false) == nothing
        @test plot(b, false, plotting = false) == nothing

    end

    @testset "copula and bivariate plots" begin

        a = πCop()

        PyPlot.ioff()

        @test plot(a) == nothing
        @test plot(a, name = "same") == nothing
        @test plot(a, fontsize = 18, alpha = 0.7) == nothing
        @test plotStep(a) == nothing

        b = N(0,1)

        j1 = bivpbox(b, b, a);

        @test plot(j1) == nothing

        samps = ProbabilityBoundsAnalysis.sample(j1, 10)

        using PyCall
        out = scatter(samps)
        @test  @isdefined out

    end

    PyPlot.ion()
end
