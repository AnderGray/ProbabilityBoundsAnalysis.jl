######
# This file is part of the ProbabilityBoundsAnalysis.jl package.
#
#   Pbox constructor and methods test
#
######


@testset "Pbox constructor" begin

    x = pbox()
    @test isvacuous(x)

    min = 1; max = 2; name = "test";
    shape = "testShape";

    n = parametersPBA.steps;

    x = pbox(min, max, name = name, shape = shape);

    @test x isa pbox
    @test isinterval(x)
    @test ispbox(x)

    @test interval(x.ml,x.mh) == interval(min,max)
    @test interval(x.vl, x.vh) == interval(0,0.25)

    @test x.name == name
    @test x.shape == shape

    @test all(x.u .== ones(n) * min)
    @test all(x.d .== ones(n) * max)

    x = pbox([0,1]);

    @test !isinterval(x)
    @test x.ml == x.mh
    @test 0.5 ∈ mean(x)
    @test x.vh == x.vl

    x = pbox(interval(0,4), ml=1, mh = 3, vl = 0.2, vh = 2)

    @test isinterval(x)
    @test x.ml == 1;
    @test x.mh == 3;
    @test x.vl == 0.2;
    @test x.vh == 2;

    y = pbox(x)
    @test y == x

    x = pbox.(0,1:3)
    @test all(isinterval.(x))

    x = makepbox(1,2,3)
    @test all(isscalar.(x))

    x = makepbox.(interval.(0,1:3))
    @test all(isinterval.(x))

    x = normal(0,1)
    @test x(0) == cdf(x, 0)


end

@testset "Pbox range" begin
    # univariate pbox
    u = uniform(interval(-1, 0), interval(0, 1))
    @test range(u) == interval(-1, 1)
end

@testset "Pbox mass" begin
    # univariate pbox

    x = interval(1.12, 1.15)
    pb = pbox(mince(x, 10))
    m = mass(pb, 1.14, 100.0)

    @test m.lo ≈ 0.295 && m.hi ≈ 0.395

    @test 0 ∈ mass(pb, interval(0, 1))
    @test 1 ∈ mass(pb, interval(0, 2))

    @test ~isempty(mass(pb, 0, 1.3) ∩ cdf(pb, 1.3))
    @test ~isempty(mass(pb, 1.3, 100) ∩ (1-cdf(pb, 1.3)))

    x1 = pbox(0, 1)
    @test mass(x1, 0.2, 0.8) == 0 .. 1
    @test mass(x1, -1, 0.8) == 0 .. 1
    @test mass(x1, 0.8, 2) == 0 .. 1
    @test 1 ∈ mass(x1, -1, 2)
    @test 0 ∈ mass(x1, -10, -2)

    @test mass(x1, 0.2, 0.8) == 0 .. 1
    @test mass(x1, -1, 0.8) == 0 .. 1
    @test mass(x1, 0.8, 2) == 0 .. 1
    @test 1 ∈ mass(x1, -1, 2)
    @test 0 ∈ mass(x1, -10, -2)

end


@testset "Pbox from interval data" begin

    ##
    # Checks if belief and plausibility from
    # DSS is the same (or larger than) the p-box
    ##
    function test_measure(Fe, masses, U)
        ids_lo = Fe .⊆ U
        ids_hi = .~isempty.(Fe .∩ U)

        P = interval(sum(masses[ids_lo]) , sum(masses[ids_hi]))
        pb = mixture(Fe, masses)

        P1_p = mass(pb, U)

        return P ⊆ P1_p
    end

    U1 = interval(0.5, 1.7)
    U2 = interval(1.5, 2)
    U3 = interval(1.7, 2)

    x1 = interval(1, 2)
    x2 = interval(1.5, 2.5)

    Fe1 = [x1, x2]

    # Check for different masses
    masses1 = [0.5, 0.5]
    masses2 = [0.3, 0.7]
    masses3 = [0.7, 0.3]
    masses4 = [0, 1]

    @test test_measure(Fe1, masses1, U1)
    @test test_measure(Fe1, masses1, U2)
    @test test_measure(Fe1, masses1, U3)

    @test test_measure(Fe1, masses2, U1)
    @test test_measure(Fe1, masses2, U2)
    @test test_measure(Fe1, masses2, U3)

    @test test_measure(Fe1, masses3, U1)
    @test test_measure(Fe1, masses3, U2)
    @test test_measure(Fe1, masses3, U3)

    @test test_measure(Fe1, masses4, U1)
    @test test_measure(Fe1, masses4, U2)
    @test test_measure(Fe1, masses4, U3)

    x3 = interval(0, 0.6)
    x4 = interval(0, 3)
    x5 = interval(2.3, 2.4)
    x6 = interval(1.5, 2)

    Fe2 = [x1, x2, x3, x4, x5, x6]

    masses1 = ones(length(Fe2)) ./length(Fe2)
    masses2 = [1, 2, 3, 4, 5, 6]
    masses2 = masses2 ./sum(masses2)

    @test test_measure(Fe2, masses1, U1)
    @test test_measure(Fe2, masses1, U2)
    @test test_measure(Fe2, masses1, U3)

    @test test_measure(Fe2, masses2, U1)
    @test test_measure(Fe2, masses2, U2)
    @test test_measure(Fe2, masses2, U3)


    pb_mix = mixture(Fe2, masses2)
    pb_ = pbox(Fe2, masses2)

    @test pb_mix.d == pb_.d
    @test pb_mix.u == pb_.u
    @test mean(pb_mix) == mean(pb_)
    @test var(pb_mix) == var(pb_)

end

@testset "p-box from 2nd order Monte Carlo" begin

    Noutter = 500;
    Ninner = 1000;

    samps = zeros(Noutter, Ninner)

    [samps[i,:] = randn(Ninner) for i = 1:Noutter]

    pb = pbox(samps)

    pb_vec = Vector{pbox}(undef, Noutter)
    for i = 1:Noutter
        edf = zeros(1, Ninner);
        edf[1,:] = samps[i,:];
        pb_vec[i] = pbox(edf)
    end

    @test all(pb_vec .⊆ pb)

end

##
# pbox functions, cdf, mean, checkmoments...
##

##
# Interpolations
##
