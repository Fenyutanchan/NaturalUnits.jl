@testset "unwrap_dimensionless_EU" begin
    @testset "dimensionless → bare value" begin
        u = eV(3.0, 0)
        r = unwrap_dimensionless_EU(u)
        @test r === 3.0
        @test r isa Float64
    end

    @testset "dimensioned → pass through unchanged" begin
        u = eV(3.0, 1)
        r = unwrap_dimensionless_EU(u)
        @test r === u
    end

    @testset "Number → identity" begin
        @test unwrap_dimensionless_EU(42) === 42
        @test unwrap_dimensionless_EU(3.14) === 3.14
        @test unwrap_dimensionless_EU(1 + 2im) === 1 + 2im
    end

    @testset "various numeric types, dimensionless" begin
        @test unwrap_dimensionless_EU(eV(7, 0)) === 7
        @test unwrap_dimensionless_EU(MeV(2.5, 0)) === 2.5
    end

    @testset "rational dimension zero is dimensionless" begin
        u = eV(4.0, 0 // 1)
        @test unwrap_dimensionless_EU(u) === 4.0
    end
end

@testset "_unwrap_dimensionless_EU" begin
    @testset "dimensionless → bare value" begin
        u = eV(2.7, 0)
        r = NaturalUnits._unwrap_dimensionless_EU(u)
        @test r === 2.7
    end

    @testset "dimensioned → ArgumentError" begin
        @test_throws ArgumentError NaturalUnits._unwrap_dimensionless_EU(eV(1.0, 1))
        @test_throws ArgumentError NaturalUnits._unwrap_dimensionless_EU(eV(1.0, 2))
    end
end

@testset "log" begin
    @testset "dimensionless" begin
        @test log(eV(1.0, 0)) === log(1.0)
        @test log(eV(ℯ, 0)) ≈ 1.0 atol = 1e-15
        @test log(eV(2.0, 0)) ≈ log(2.0)
    end

    @testset "dimensioned → ArgumentError" begin
        @test_throws ArgumentError log(eV(1.0, 1))
    end
end

@testset "log2" begin
    @testset "dimensionless" begin
        @test log2(eV(1.0, 0)) === log2(1.0)
        @test log2(eV(8.0, 0)) ≈ 3.0
    end

    @testset "dimensioned → ArgumentError" begin
        @test_throws ArgumentError log2(eV(1.0, 1))
    end
end

@testset "log10" begin
    @testset "dimensionless" begin
        @test log10(eV(1.0, 0)) === log10(1.0)
        @test log10(eV(100.0, 0)) ≈ 2.0
    end

    @testset "dimensioned → ArgumentError" begin
        @test_throws ArgumentError log10(eV(1.0, 1))
    end
end

@testset "exp" begin
    @testset "dimensionless" begin
        @test exp(eV(0.0, 0)) === 1.0
        @test exp(eV(1.0, 0)) ≈ ℯ atol = 1e-15
        @test exp(eV(2.0, 0)) ≈ exp(2.0)
    end

    @testset "dimensioned → ArgumentError" begin
        @test_throws ArgumentError exp(eV(0.0, 1))
    end
end

@testset "exp2" begin
    @testset "dimensionless" begin
        @test exp2(eV(0.0, 0)) === 1.0
        @test exp2(eV(3.0, 0)) ≈ 8.0
    end

    @testset "dimensioned → ArgumentError" begin
        @test_throws ArgumentError exp2(eV(0.0, 1))
    end
end

@testset "exp10" begin
    @testset "dimensionless" begin
        @test exp10(eV(0.0, 0)) === 1.0
        @test exp10(eV(2.0, 0)) ≈ 100.0
    end

    @testset "dimensioned → ArgumentError" begin
        @test_throws ArgumentError exp10(eV(0.0, 1))
    end
end

@testset "exp/log round-trip" begin
    u = eV(1.5, 0)
    @test exp(log(u)) ≈ 1.5
    @test log(exp(u)) ≈ 1.5

    @test exp2(log2(u)) ≈ 1.5
    @test exp10(log10(u)) ≈ 1.5
end

@testset "log/exp with other prefix types" begin
    @test log(MeV(1.0, 0)) === log(1.0)
    @test exp(GeV(0.0, 0)) === 1.0
    @test log10(keV(100.0, 0)) ≈ 2.0
end

@testset "log/exp with integer-valued dimensionless" begin
    @test log(eV(1, 0)) ≈ 0.0
    @test exp(eV(0, 0)) ≈ 1.0
end
