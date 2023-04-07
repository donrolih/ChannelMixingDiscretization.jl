using ChannelMixingDiscretization
using Test

@testset "Star model mapping" begin
    # Write your tests here.
    S = Float64
    J = 35
    Nz = 10
    n = 2
    E = zeros(S, J, Nz, n, n)
    T = zeros(S, J, Nz, n, n)
    zs = collect(range(0.5/Nz, 1-0.5/Nz, Nz))
    @test ChannelMixingDiscretization.twistingparameters(Nz) == zs
end