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
