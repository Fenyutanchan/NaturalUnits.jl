@testset "eV constructors" begin
    @testset "default constructor" begin
        u = eV()
        @test EUval(u) == 1
        @test EUdim(u) == 1
    end

    @testset "single-argument constructor" begin
        u = eV(3.0)
        @test EUval(u) == 3.0
        @test EUdim(u) == 1

        u_int = eV(5)
        @test EUval(u_int) == 5
        @test EUdim(u_int) == 1
    end

    @testset "two-argument constructor" begin
        u = eV(2.0, 2)
        @test EUval(u) == 2.0
        @test EUdim(u) == 2

        u_rat = eV(3.0, 1 // 2)
        @test EUval(u_rat) == 3.0
        @test EUdim(u_rat) == 1 // 2
    end

    @testset "dimension zero returns bare value" begin
        val = eV(42.0, 0)
        @test val isa eV
        @test EUval(val) === 42.0
        @test EUdim(val) == 0

        val_int = eV(7, 0)
        @test val_int isa eV
        @test EUval(val_int) === 7
        @test EUdim(val_int) == 0
    end
end

@testset "prefix type constructors (keV, MeV, GeV, TeV)" begin
    for (T, sym) in [(keV, :keV), (MeV, :MeV), (GeV, :GeV), (TeV, :TeV)]
        @testset "$sym" begin
            u_default = T()
            @test EUval(u_default) == 1
            @test EUdim(u_default) == 1

            u_val = T(2.5)
            @test EUval(u_val) == 2.5
            @test EUdim(u_val) == 1

            u_dim = T(3.0, 2)
            @test EUval(u_dim) == 3.0
            @test EUdim(u_dim) == 2

            val_zero = T(10.0, 0)
            @test val_zero isa T
            @test EUval(val_zero) === 10.0
            @test EUdim(val_zero) == 0
        end
    end
end

@testset "EUdim and EUval accessors" begin
    @testset "EUdim" begin
        @test EUdim(eV(1, 3)) == 3
        @test EUdim(eV(1, 1 // 2)) == 1 // 2
        @test EUdim(42) == 0
        @test EUdim(3.14) == 0
    end

    @testset "EUval" begin
        @test EUval(eV(7.0, 1)) == 7.0
        @test EUval(42) == 42
        @test EUval(3.14) == 3.14
    end

    @testset "EUval with type" begin
        u = MeV(2.0, 1)
        val_in_eV = EUval(eV, u)
        @test val_in_eV ≈ 2.0e6
    end
end

@testset "one, oneunit and zero" begin
    @testset "one is dimensionless multiplicative identity" begin
        u = eV(5.0, 3)
        @test one(u) isa eV
        @test EUval(one(u)) == 1
        @test EUdim(one(u)) == 0
        @test one(u) * u == u
        @test u * one(u) == u

        @test one(eV) isa eV
        @test EUval(one(eV)) == 1
        @test EUdim(one(eV)) == 0

        @test one(MeV) isa MeV
        @test EUval(one(MeV)) == 1
        @test EUdim(one(MeV)) == 0

        @test one(GeV) isa GeV
        @test EUval(one(GeV)) == 1
        @test EUdim(one(GeV)) == 0

        # one(T) and one(u::T) must be equal as multiplicative identities
        @test one(u) == one(eV)
    end

    @testset "oneunit preserves dimension" begin
        u = eV(5.0, 3)
        @test oneunit(u) isa eV
        @test EUval(oneunit(u)) == 1
        @test EUdim(oneunit(u)) == 3

        @test oneunit(eV) isa eV
        @test EUval(oneunit(eV)) == 1
        @test EUdim(oneunit(eV)) == 1

        @test oneunit(MeV) isa MeV
        @test EUval(oneunit(MeV)) == 1
        @test EUdim(oneunit(MeV)) == 1
    end

    @testset "zero" begin
        u = eV(5.0, 3)
        @test zero(u) isa eV
        @test EUval(zero(u)) == 0
        @test EUdim(zero(u)) == 3
        @test zero(u) + u == u
        @test u + zero(u) == u

        # zero preserves dimension across instances
        v = GeV(2.0, 4)
        @test zero(v) isa GeV
        @test EUdim(zero(v)) == 4

        # zero on the bare type is intentionally undefined: there is no
        # canonical additive identity without a known dimension
        @test_throws MethodError zero(eV)
        @test_throws MethodError zero(MeV)
        @test_throws MethodError zero(GeV)
        @test_throws MethodError zero(EnergyUnit)
    end

    @testset "numeric zero is additive identity at any dimension" begin
        u = eV(5.0, 3)
        @test u + 0 == u
        @test 0 + u == u
        @test u - 0 == u
        @test 0 - u == -u

        v = GeV(7.0, 2)
        @test v + 0 === v
        @test v + 0.0 == v
        @test v + 0 // 1 == v

        # non-zero numbers must still respect dimensions
        @test_throws ArgumentError u + 1
        @test_throws ArgumentError 1 + u
        @test_throws ArgumentError u - 1
    end
end
