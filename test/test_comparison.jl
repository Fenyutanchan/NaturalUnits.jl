@testset "equality ==" begin
    @testset "same value and dimension" begin
        @test eV(3.0, 1) == eV(3.0, 1)
    end

    @testset "different value" begin
        @test eV(3.0, 1) != eV(4.0, 1)
    end

    @testset "different dimension" begin
        @test eV(3.0, 1) != eV(3.0, 2)
    end

    @testset "cross-type equality" begin
        @test eV(1.0e6, 1) == MeV(1.0, 1)
    end

    @testset "cross-type inequality" begin
        @test eV(1.0, 1) != MeV(1.0, 1)
    end
end

@testset "isless" begin
    @testset "same dimension" begin
        @test isless(eV(1.0, 1), eV(2.0, 1))
        @test !isless(eV(2.0, 1), eV(1.0, 1))
        @test !isless(eV(1.0, 1), eV(1.0, 1))
    end

    @testset "different dimension throws" begin
        @test_throws ArgumentError isless(eV(1.0, 1), eV(1.0, 2))
    end

    @testset "cross-type comparison" begin
        @test isless(MeV(1.0, 1), GeV(1.0, 1))
        @test isless(keV(1.0, 1), MeV(1.0, 1))
    end
end

@testset "predicates" begin
    @testset "isinf" begin
        @test isinf(eV(Inf, 1))
        @test !isinf(eV(1.0, 1))
    end

    @testset "isnan" begin
        @test isnan(eV(NaN, 1))
        @test !isnan(eV(1.0, 1))
    end

    @testset "iszero" begin
        @test iszero(eV(0.0, 1))
        @test !iszero(eV(1.0, 1))
    end
end

@testset "== vs isless asymmetry (intentional)" begin
    a = eV(1.0, 1)
    b = eV(1.0, 2)
    @test (a == b) === false
    @test_throws ArgumentError isless(a, b)
end

@testset "isequal & hash" begin
    a = eV(1.0e6, 1)
    b = MeV(1.0, 1)
    @test a == b
    @test isequal(a, b)
    @test hash(a) == hash(b)

    c = eV(1.0, 1)
    d = eV(1.0, 2)
    @test !isequal(c, d)
    @test hash(c) != hash(d)

    s = Set([a, b, c])
    @test length(s) == 2

    dct = Dict(a => "x")
    @test dct[b] == "x"
end

@testset "NaN equality semantics" begin
    n1 = eV(NaN, 1)
    n2 = eV(NaN, 1)
    @test (n1 == n2) === false
    @test isequal(n1, n2)
    @test hash(n1) == hash(n2)
end
