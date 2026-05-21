const nu_eV = NaturalUnit(eV)
const nu_GeV = NaturalUnit(GeV)

@testset "1 Joule in natural units" begin
    J_in_eV = nu_eV.J
    @test EUdim(J_in_eV) == 1
    @test EUval(J_in_eV) ≈ 1.0 / 1.602176634e-19 rtol = 1e-6

    J_in_GeV = nu_GeV.J
    @test EUdim(J_in_GeV) == 1
    @test EUval(J_in_GeV) ≈ 1.0 / 1.602176634e-10 rtol = 1e-6
end

@testset "1 meter in natural units" begin
    m_in_GeV = nu_GeV.m
    @test EUdim(m_in_GeV) == -1
    @test EUval(m_in_GeV) ≈ 5.0677307e15 rtol = 1e-4
end

@testset "1 centimeter in natural units" begin
    cm_in_GeV = nu_GeV.cm
    @test EUdim(cm_in_GeV) == -1
    @test EUval(cm_in_GeV) ≈ EUval(nu_GeV.m) * 1e-2 rtol = 1e-14
end

@testset "1 second in natural units" begin
    s_in_GeV = nu_GeV.s
    @test EUdim(s_in_GeV) == -1
    @test EUval(s_in_GeV) ≈ 1.5192674e24 rtol = 1e-4
end

@testset "1 kilogram in natural units" begin
    kg_in_GeV = nu_GeV.kg
    @test EUdim(kg_in_GeV) == 1
    @test EUval(kg_in_GeV) ≈ 5.6095886e26 rtol = 1e-4
end

@testset "1 gram in natural units" begin
    g_in_GeV = nu_GeV.g
    @test EUdim(g_in_GeV) == 1
    @test EUval(g_in_GeV) ≈ EUval(nu_GeV.kg) * 1e-3 rtol = 1e-14
end

@testset "1 Kelvin in natural units" begin
    K_in_GeV = nu_GeV.K
    @test EUdim(K_in_GeV) == 1
    @test EUval(K_in_GeV) ≈ 8.6173333e-14 rtol = 1e-4
end

@testset "cross-check: ℏc = 197.327 MeV·fm" begin
    m_in_MeV = NaturalUnit(MeV).m
    hbar_c_MeV_fm = 1.0 / EUval(m_in_MeV) * 1e15
    @test hbar_c_MeV_fm ≈ 197.3269804 rtol = 1e-6
end

@testset "cross-check: 1 GeV/c² in kg" begin
    GeV_in_kg = 1.78266192e-27
    inv_kg_in_GeV = 1.0 / GeV_in_kg
    @test EUval(nu_GeV.kg) ≈ inv_kg_in_GeV rtol = 1e-4
end
