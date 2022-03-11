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



        pl1 = plot(a, plotting = false)
        pl2 = plot(a, name = "same", plotting = false)
        pl3 = plot(a, col = "red", plotting = false)
        pl4 = plot(a, alpha = 0.5, plotting = false)
        pl5 = plot(a, fontsize = 20, plotting = false)
        pl6 = plot(a, false, plotting = false)

        @test @isdefined pl1
        @test @isdefined pl2
        @test @isdefined pl3
        @test @isdefined pl4
        @test @isdefined pl5
        @test @isdefined pl6

        b = interval(2, 3)

        pl1_i = plot(b, plotting = false)
        pl2_i = plot(b, name = "same", plotting = false)
        pl3_i = plot(b, col = "red", plotting = false)
        pl4_i = plot(b, alpha = 0.5, plotting = false)
        pl5_i = plot(b, fontsize = 20, plotting = false)
        pl6_i = plot(b, false, plotting = false)

        @test @isdefined pl1_i
        @test @isdefined pl2_i
        @test @isdefined pl3_i
        @test @isdefined pl4_i
        @test @isdefined pl5_i
        @test @isdefined pl6_i

    end

    @testset "copula and bivariate plots" begin

        a = πCop()

        PyPlot.ioff()

        pl1 = plot(a)
        pl2 = plot(a, name = "same")
        pl3 = plot(a, fontsize = 18, alpha = 0.7)
        pl4 = plotStep(a)

        @test @isdefined pl1
        @test @isdefined pl2
        @test @isdefined pl3
        @test @isdefined pl4

        b = N(0,1)

        j1 = bivpbox(b, b, a);

        plj_1 = plot(j1)

        @test @isdefined plj_1

        samps = ProbabilityBoundsAnalysis.sample(j1, 10)

        plj_2 = scatter(samps)

        @test @isdefined plj_2
    end

    PyPlot.ion()
end
