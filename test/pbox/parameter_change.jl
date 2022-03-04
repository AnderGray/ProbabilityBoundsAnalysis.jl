@testset "Change parameters" begin

    Nsteps = 100
    ProbabilityBoundsAnalysis.setSteps(Nsteps)

    x = U(0, 1)

    @test ProbabilityBoundsAnalysis.parametersPBA.steps == Nsteps
    @test length(x.d) == Nsteps
    @test x.n == Nsteps

    newBot = 0.1
    ProbabilityBoundsAnalysis.setBOt(newBot)
    @test ProbabilityBoundsAnalysis.parametersPBA.bOt == newBot
    @test ProbabilityBoundsAnalysis.iii()[1] == newBot

    newTop = 0.1
    ProbabilityBoundsAnalysis.setTOp(newTop)
    @test ProbabilityBoundsAnalysis.parametersPBA.tOp == newTop
    @test ProbabilityBoundsAnalysis.jjj()[end] == newTop

    newVerbose = 1;
    ProbabilityBoundsAnalysis.setVerbose(newVerbose)
    @test ProbabilityBoundsAnalysis.parametersPBA.verbose== newVerbose

end
