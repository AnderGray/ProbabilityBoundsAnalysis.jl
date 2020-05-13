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

    n = ProbabilityBoundsAnalysis.steps;

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
    @test x.ml == 0.5
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

    
end

##
# pbox functions, cdf, mean, checkmoments...
##

##
# Interpolations
##

