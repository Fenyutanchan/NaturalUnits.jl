@testset "broadcastable" begin
    u = eV(3.0, 1)
    b = Broadcast.broadcastable(u)
    @test b isa Base.RefValue{<:EnergyUnit}
    @test b[] === u
end

@testset "broadcasting" begin
    u = eV(3.0, 1)

    # Scalar broadcasting over a scalar argument unwraps to the plain value.
    @test EUval.(u) === 3.0

    # Broadcasting with an array: the unit acts as a scalar and is reused.
    arr = [1.0, 2.0, 3.0]
    @test (arr .* u) == [eV(3.0, 1), eV(6.0, 1), eV(9.0, 1)]
end
