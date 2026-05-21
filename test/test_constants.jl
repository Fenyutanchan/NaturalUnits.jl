const nu_for_constants = NaturalUnit(GeV)

@testset "Newton's gravitational constant G" begin
    G = nu_for_constants.G_N
    @test EUdim(G) == -2
end

@testset "reduced Planck mass M_Pl" begin
    M_Pl = nu_for_constants.M_Pl
    @test EUdim(M_Pl) == 1
    @test EUval(M_Pl) ≈ 2.4353e18 rtol = 1e-2
end

@testset "Planck mass m_Pl" begin
    m_Pl = nu_for_constants.m_Pl
    @test EUdim(m_Pl) == 1
    @test EUval(m_Pl) ≈ 1.2209e19 rtol = 1e-2
end

@testset "self-consistency: M_Pl = 1/√(8πG)" begin
    G = nu_for_constants.G_N
    M_Pl = nu_for_constants.M_Pl
    @test EUval(M_Pl) ≈ 1.0 / sqrt(8 * π * EUval(G)) rtol = 1e-10
end

@testset "self-consistency: m_Pl = 1/√G" begin
    G = nu_for_constants.G_N
    m_Pl = nu_for_constants.m_Pl
    @test EUval(m_Pl) ≈ 1.0 / sqrt(EUval(G)) rtol = 1e-10
end

@testset "self-consistency: M_Pl < m_Pl" begin
    @test EUval(nu_for_constants.M_Pl) < EUval(nu_for_constants.m_Pl)
end

@testset "self-consistency: G dimensionally correct" begin
    m = nu_for_constants.m
    kg = nu_for_constants.kg
    s = nu_for_constants.s
    G_dimension = EUdim(m^3) - EUdim(kg) - EUdim(s^2)
    @test G_dimension == EUdim(nu_for_constants.G_N)
end

@testset "self-consistency: 1 J = 1 kg·m²/s²" begin
    J = nu_for_constants.J
    kg = nu_for_constants.kg
    m = nu_for_constants.m
    s = nu_for_constants.s
    @test EUval(J) ≈ EUval(kg) * EUval(m)^2 / EUval(s)^2 rtol = 1e-8
end
