@testset "convert between prefix types" begin
    @testset "MeV -> eV" begin
        u = MeV(1.0, 1)
        c = convert(eV, u)
        @test c isa eV
        @test EUval(c) == 1.0e6
        @test EUdim(c) == 1
    end

    @testset "GeV -> eV" begin
        u = GeV(1.0, 1)
        c = convert(eV, u)
        @test EUval(c) == 1.0e9
        @test EUdim(c) == 1
    end

    @testset "keV -> eV" begin
        u = keV(2.0, 1)
        c = convert(eV, u)
        @test EUval(c) == 2.0e3
        @test EUdim(c) == 1
    end

    @testset "TeV -> eV" begin
        u = TeV(1.0, 1)
        c = convert(eV, u)
        @test EUval(c) == 1.0e12
        @test EUdim(c) == 1
    end

    @testset "eV -> MeV" begin
        u = eV(1.0e6, 1)
        c = convert(MeV, u)
        @test c isa MeV
        @test EUval(c) == 1.0
        @test EUdim(c) == 1
    end

    @testset "eV -> GeV" begin
        u = eV(3.0e9, 1)
        c = convert(GeV, u)
        @test EUval(c) == 3.0
        @test EUdim(c) == 1

        u2 = eV(3.0e18, 2)
        c2 = convert(GeV, u2)
        @test EUval(c2) == 3.0
        @test EUdim(c2) == 2
    end

    @testset "GeV -> MeV (via eV)" begin
        u = GeV(1.0, 1)
        c = convert(MeV, u)
        @test c isa MeV
        @test EUval(c) == 1.0e3
        @test EUdim(c) == 1
    end

    @testset "dimension preserved in conversion" begin
        for dim in [1, 2, -1, 1 // 2, 3 // 2]
            u = GeV(2.0, dim)
            c = convert(eV, u)
            @test EUdim(c) == dim
        end
    end
end

@testset "convert same type is identity" begin
    u = eV(5.0, 2)
    @test convert(eV, u) === u

    u2 = MeV(3.0, 1)
    @test convert(MeV, u2) === u2
end

@testset "convert Number → EnergyUnit" begin
    c = convert(eV, 42.0)
    @test c isa eV
    @test EUval(c) === 42.0
    @test EUdim(c) == 0

    c2 = convert(EnergyUnit, 42.0)
    @test c2 isa eV
    @test EUval(c2) === 42.0
    @test EUdim(c2) == 0

    c3 = convert(GeV, 7)
    @test c3 isa GeV
    @test EUval(c3) === 7
    @test EUdim(c3) == 0
end

@testset "promote_rule" begin
    @test promote_rule(MeV, eV) == eV
    @test promote_rule(GeV, keV) == eV
    @test promote_rule(TeV, MeV) == eV
    @test promote_rule(eV, eV) == eV
end

@testset "mixed-type operations use promote" begin
    r = keV(1) + MeV(1)
    @test r isa eV
    @test EUval(r) == 1.0e3 + 1.0e6
    @test EUdim(r) == 1

    r2 = GeV(1) * keV(1)
    @test r2 isa eV
    @test EUval(r2) == 1.0e9 * 1.0e3
    @test EUdim(r2) == 2
end

@testset "convert_EnergyUnit_value_type" begin
    u = eV(3, 2)
    u_float = convert_EnergyUnit_value_type(Float64, u)
    @test EUval(u_float) isa Float64
    @test EUval(u_float) == 3.0
    @test EUdim(u_float) == 2
end
