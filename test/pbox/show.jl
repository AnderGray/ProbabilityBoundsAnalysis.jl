######
# This file is part of the ProbabilityBoundsAnalysis.jl package.
#
#   Show tests
#
######

replstr(x, kv::Pair...) = sprint((io,x) -> show(IOContext(io, :limit => true, :displaysize => (24, 80), kv...), MIME("text/plain"), x), x)
showstr(x, kv::Pair...) = sprint((io,x) -> show(IOContext(io, :limit => true, :displaysize => (24, 80), kv...), x), x)

@testset "pbox show" begin

    x = U(0,1)
    @test replstr(x) == "Pbox: \t  ~ uniform ( range=[0.0, 1.0], mean=0.5, var=0.083333)"

    x2 = U(0, interval(1, 2))
    @test replstr(x2) == "Pbox: \t  ~ uniform ( range=[0.0, 2.0], mean=[0.5, 1.0], var=[0.083333, 0.33333])"

end

@testset "copula show" begin

    x = πCop()
    @test replstr(x) == "Copula ~ π()"

    x2 = GauCopula()
    @test replstr(x2) == "Copula ~ Gaussian(r = 0)"

    x3 = GauCopula(interval(0, 1))
    @test replstr(x2) == "Imp Copula ~ Gaussian(r = [0, 1])"

end

@testset "bivariate p-box show" begin

    x1 = U(0, 1)
    x2 = U(0, 1)

    c = πCop()
    j1 = c(x1, x2)

    @test replstr(j1) == "BivPbox ~ π( uniform(mean = 0.5, var = 0.08333333333333333), uniform(mean = 0.5, var = 0.08333333333333333))"

    y1 = U(0..1, 1..2)
    y2 = U(0..1, 1..2)
    c2 = Frechet()

    j2 = c2(y1, y2)
    @test replstr(j2) == "BivPbox ~ Arb( uniform(mean = [0.5, 1.5], var = [0, 0.333334]), uniform(mean = [0.5, 1.5], var = [0, 0.333334]))"

end
