######
# This file is part of the ProbabilityBoundsAnalysis.jl package.
#
#   Tests special functions
#
######

@testset "special functions" begin

    @test interval(interval(0,1), interval(1,2)) == interval(0, 2)

    @test left(1) == 1
    @test right(1) == 1

    int = interval(0, 1);
    pb = U(0,1)

    @test ispbox(pb)

    @test left(int) == 0
    @test right(int) == 1

    @test left(pb) == 0
    @test right(pb) == 1

    @test all(lefts(pb) .== pb.u)
    @test all(rights(pb) .== pb.d)

    @test isinterval(int)
    @test isinterval(U(0,0))

    @test isscalar(0)
    @test isscalar(interval(0,0))
    @test isscalar(U(0,0))

    @test isvacuous(interval(-Inf, Inf))
    @test isvacuous(pbox(-Inf, Inf))

    @test straddles(interval(-1, 1))
    @test straddles(U(-1, 1))

    @test straddles(interval(0, 1))
    @test straddles(U(0, 1))

    @test straddles(interval(-1, 0))
    @test straddles(U(-1, 0))

    @test ~straddles(interval(2, 3))
    @test ~straddles(U(2, 3))

    @test ~straddles(interval(-3, -2))
    @test ~straddles(U(-3, -2))

    @test straddlingzero(interval(-1, 1))
    @test straddlingzero(U(-1, 1))

    @test ~straddlingzero(interval(-1, 0))
    @test ~straddlingzero(U(-1, 0))
    @test ~straddlingzero(interval(0, 1))
    @test ~straddlingzero(U(0, 1))

    @test issubset([0,0], interval(-1, 1) × interval(-1, 1))
    @test issubset( -1 .. 1, interval(-1, 1) × interval(-1, 1))
    @test issubset(interval(-1, 1) × interval(-1, 1),  -1 .. 1)

    @test intersect(2.0, interval(-1,3)) == 2
    @test 2 ∩ interval(-1,3) == 2
    @test isempty( 4 ∩ interval(-1,3))
    @test isempty(interval(-1,3) ∩ 4.0)


    Ints = mince(0.. 1, 40)
    @test no_nesting(Ints)
    @test ~no_nesting([Ints; 0.2..0.3])

    @test touching(1, 1.0)
    @test touching(1, interval(0.5,1))
    @test touching(interval(0,0.5), interval(0.5,1))
    @test touching(1, U(0.5,1))
    @test touching(U(0,0.5), U(0.5,1))
    @test touching(interval(0,0.5), U(0.5,1))

    @test ~touching(2, 3)
    @test ~touching(interval(0,0.5), 3)
    @test ~touching(interval(0,0.5), 2..3)
    @test ~touching(interval(0,0.5), U(2, 3))


end
