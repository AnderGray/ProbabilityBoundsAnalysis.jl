######
# This file is part of the ProbabilityBoundsAnalysis.jl package.
#
#   Comparisons tests
#
######

@testset "Comparisons" begin

    @testset "scalar comparisons" begin

        a = U(0,1)


        @test a <= 1
        @test a >= 0
        @test a < 1.1
        @test a > -0.1

        @test 1 >= a
        @test 0 <= a
        @test 1.1 > a
        @test -0.1 < a

        @test 0.6 ∈ (a < 0.6)
        @test 0.4 ∈ (a > 0.6)

        @test 0.6 ∈ (a <= 0.6)
        @test 0.4 ∈ (a >= 0.6)

        @test 0.6 ∈ (0.6 > a)
        @test 0.4 ∈ (0.6 < a)

        @test 0.6 ∈ (0.6 >= a)
        @test 0.4 ∈ (0.6 <= a)

    end

    @testset "interval comparisons" begin

        a = U(0,1)

        @test a < interval(1.001,2)
        @test a <= interval(1.001,2)
        @test a >= interval(-2, -0.001)
        @test a > interval(-2, -0.001)

        @test interval(1.001,2) > a
        @test interval(1.001,2) >= a
        @test interval(-2, -0.001) <= a
        @test interval(-2, -0.001) < a

        b = interval(0.6, 1)

        @test b ⊆ ( a < b)
        @test b ⊆ ( a <= b)
        @test 1 - b ⊆ ( a > b)
        @test 1 - b ⊆ ( a >= b)

        @test b ⊆ ( b > a)
        @test b ⊆ ( b >= a)
        @test 1 - b ⊆ ( b < a)
        @test 1 - b ⊆ ( b <= a)

    end

    @testset "p-box comparisons" begin

        a1 = U(0, 1);
        a2 = U(1, 2);
        a3 = U(2, 3);
        a4 = U(0.5, 1.5);

        @test a1 <= a2
        @test 1 ∈ (a1 < a2)
        @test 0 ∈ (a2 <= a1)

        @test ~(a1 >= a3)
        @test ~(a1 > a3)

        @test interval(0.5, 1) ⊆ (a1 < a4)
        @test interval(0.5, 1) ⊆ (a4 > a1)

        @test interval(0, 0.5) ⊆ (a1 > a4)
        @test interval(0, 0.5) ⊆ (a4 < a1)

        @test <(a1, a4, corr = 1)
        @test >(a4, a1, corr = 1)

        @test <=(a1, a4, corr = 1)
        @test >=(a4, a1, corr = 1)

        @test 0.75 ∈ <(a1, a4, corr = -1)
        @test 0.75 ∈ >(a4, a1, corr = -1)

        @test 0.875 ∈ <(a1, a4, corr = 0)
        @test 0.875 ∈ >(a4, a1, corr = 0)


    end

    @testset "Subset checks" begin

        x1 = U(0, 2)
        x2 = U(-1..0, 1..2)

        @test x1 ⊆ x2
        @test ~(x1 ⊂ x2)

        x3 = U(-0.5, 1.5)
        @test x3 ⊂ x2

        x4 = U(-0.9 .. -0.5, 1.1..1.7)

        @test x4 ⊂ x2
        @test x4 ⊆ x2

        x5 = U(-1, 3)

        @test ~(x5 ⊂ x2)
        @test ~(x5 ⊆ x2)

    end

    @testset "equality checks" begin

        x1 = U(0, 2)
        x2 = U(0, 2)

        x3 = deepcopy(x2)
        x3.d[end] = Inf

        x4 = deepcopy(x2)
        x4.u[1] = -Inf

        x5 = deepcopy(x2)
        x5.shape = ""

        x6 = deepcopy(x2)
        x6.ml = 0;

        x7 = deepcopy(x2)
        x7.vl = 0

        x8 = N(0, 1)

        @test x1 == x2
        @test ~(x3 == x2)
        @test ~(x4 == x2)
        @test ~(x5 == x2)
        @test ~(x6 == x2)
        @test ~(x7 == x2)
        @test ~(x8 == x2)

    end

end
