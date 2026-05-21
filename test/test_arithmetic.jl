@testset "addition and subtraction" begin
    @testset "same type, same dimension" begin
        a = eV(2.0, 1)
        b = eV(3.0, 1)
        @test EUval(a + b) ≈ 5.0
        @test EUdim(a + b) == 1
        @test EUval(a - b) ≈ -1.0
        @test EUdim(a - b) == 1
    end

    @testset "same type, different dimension throws" begin
        a = eV(1.0, 1)
        b = eV(1.0, 2)
        @test_throws ArgumentError a + b
        @test_throws ArgumentError a - b
    end

    @testset "mixed types promote to eV" begin
        a = keV(1.0, 1)
        b = MeV(1.0, 1)
        r = a + b
        @test r isa eV
        @test EUval(r) ≈ 1.0e3 + 1.0e6
    end

    @testset "mixed types different dimension throws" begin
        a = keV(1.0, 1)
        b = MeV(1.0, 2)
        @test_throws ArgumentError a + b
    end

    @testset "unary minus" begin
        a = eV(3.0, 2)
        @test EUval(-a) == -3.0
        @test EUdim(-a) == 2
    end
end

@testset "multiplication" begin
    @testset "EnergyUnit * EnergyUnit" begin
        a = eV(2.0, 1)
        b = eV(3.0, 2)
        r = a * b
        @test EUval(r) ≈ 6.0
        @test EUdim(r) == 3
    end

    @testset "num * EnergyUnit" begin
        a = 2.0 * eV(3.0, 1)
        @test EUval(a) ≈ 6.0
        @test EUdim(a) == 1
    end

    @testset "EnergyUnit * num" begin
        a = eV(3.0, 1) * 2.0
        @test EUval(a) ≈ 6.0
        @test EUdim(a) == 1
    end

    @testset "dimension addition" begin
        a = eV(2.0, -1)
        b = eV(3.0, 2)
        @test EUdim(a * b) == 1
    end
end

@testset "division" begin
    @testset "EnergyUnit / EnergyUnit" begin
        a = eV(6.0, 3)
        b = eV(2.0, 1)
        r = a / b
        @test EUval(r) ≈ 3.0
        @test EUdim(r) == 2
    end

    @testset "num / EnergyUnit" begin
        r = 6.0 / eV(2.0, 1)
        @test EUval(r) ≈ 3.0
        @test EUdim(r) == -1
    end

    @testset "EnergyUnit / num" begin
        a = eV(6.0, 2) / 3.0
        @test EUval(a) ≈ 2.0
        @test EUdim(a) == 2
    end

    @testset "dimension subtraction" begin
        a = eV(1.0, 1)
        b = eV(1.0, 1)
        r = a / b
        @test r isa eV
        @test EUdim(r) == 0
        @test EUval(r) ≈ 1.0
    end
end

@testset "rational division //" begin
    a = eV(6, 3)
    b = eV(2, 1)
    r = a // b
    @test EUval(r) == 6 // 2
    @test EUdim(r) == 2

    r2 = 6 // eV(2, 1)
    @test EUdim(r2) == -1

    r3 = eV(6, 2) // 3
    @test EUdim(r3) == 2
end

@testset "power" begin
    @testset "integer power" begin
        a = eV(2.0, 1)
        r = a^3
        @test EUval(r) ≈ 8.0
        @test EUdim(r) == 3
    end

    @testset "rational power" begin
        a = eV(4.0, 2)
        r = a^(1 // 2)
        @test EUval(r) ≈ 2.0
        @test EUdim(r) == 1
    end

    @testset "negative power" begin
        a = eV(2.0, 2)
        r = a^-1
        @test EUval(r) ≈ 0.5
        @test EUdim(r) == -2
    end
end

@testset "inv" begin
    a = eV(4.0, 3)
    r = inv(a)
    @test EUval(r) ≈ 0.25
    @test EUdim(r) == -3
end
