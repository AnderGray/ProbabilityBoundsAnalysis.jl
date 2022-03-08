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

end
