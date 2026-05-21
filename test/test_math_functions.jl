@testset "abs" begin
    a = eV(-3.0, 2)
    r = abs(a)
    @test EUval(r) ≈ 3.0
    @test EUdim(r) == 2
end

@testset "abs2" begin
    a = eV(3.0, 1)
    r = abs2(a)
    @test EUval(r) ≈ 9.0
    @test EUdim(r) == 2
end

@testset "sqrt" begin
    a = eV(9.0, 4)
    r = sqrt(a)
    @test EUval(r) ≈ 3.0
    @test EUdim(r) == 2
end

@testset "cbrt" begin
    a = eV(27.0, 6)
    r = cbrt(a)
    @test EUval(r) ≈ 3.0
    @test EUdim(r) == 2
end

@testset "complex-valued operations" begin
    a = eV(3.0 + 4.0im, 1)

    @testset "real" begin
        r = real(a)
        @test EUval(r) == 3.0
        @test EUdim(r) == 1
    end

    @testset "imag" begin
        r = imag(a)
        @test EUval(r) == 4.0
        @test EUdim(r) == 1
    end

    @testset "conj" begin
        r = conj(a)
        @test EUval(r) == 3.0 - 4.0im
        @test EUdim(r) == 1
    end

    @testset "angle" begin
        r = angle(a)
        @test r ≈ atan(4.0, 3.0)
        @test r isa Number
    end
end
