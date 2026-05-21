@testset "zero-value operations" begin
    z = eV(0.0, 1)
    @test iszero(z)
    @test EUval(z + z) == 0.0
    @test EUval(z * eV(3.0, 2)) == 0.0
end

@testset "Inf and NaN propagation" begin
    u_inf = eV(Inf, 1)
    @test isinf(u_inf)

    u_nan = eV(NaN, 1)
    @test isnan(u_nan)

    r = eV(1.0, 1) + u_inf
    @test isinf(r)

    r2 = eV(1.0, 1) + u_nan
    @test isnan(r2)
end

@testset "dimension mismatch errors" begin
    a = eV(1.0, 1)
    b = eV(1.0, 2)

    @test_throws ArgumentError a + b
    @test_throws ArgumentError a - b
    @test_throws ArgumentError isless(a, b)
end

@testset "add_property_function" begin
    called = Ref(false)
    function __test_prop(nu::NaturalUnit)
        called[] = true
        return 42
    end

    @testset "registering a fresh property does not warn" begin
        @test_logs add_property_function(:test_prop, __test_prop)
    end

    nu = NaturalUnit(GeV)
    @test nu.test_prop == 42
    @test called[]

    @testset "overwriting an existing property warns" begin
        called2 = Ref(false)
        function __test_prop_v2(nu::NaturalUnit)
            called2[] = true
            return 99
        end
        @test_logs (:warn, r"Overwriting existing property function for test_prop"i) add_property_function(:test_prop, __test_prop_v2)
        # the new function replaces the old one
        @test nu.test_prop == 99
        @test called2[]
    end

    @testset "fallback to getfield for unregistered names" begin
        # `unit` is an actual struct field, not in the dict
        @test nu.unit === GeV
        # truly unknown names propagate the getfield error
        @test_throws FieldError nu.definitely_not_a_property
    end

    @testset "propertynames is consistent with getproperty" begin
        pn = propertynames(nu)
        # the real struct field is exposed
        @test :unit in pn
        # built-in registered properties
        for name in (:J, :m, :cm, :s, :kg, :g, :K, :G_N, :M_Pl, :m_Pl)
            @test name in pn
        end
        # user-registered property added earlier in this testset
        @test :test_prop in pn
        # every reported name is actually accessible without throwing
        for name in pn
            @test getproperty(nu, name) !== nothing
        end
        # private=true returns the same set (NaturalUnit has no private fields beyond :unit)
        @test Set(propertynames(nu, true)) == Set(pn)
    end
end

@testset "assertion macros" begin
    @testset "@check_EU_dimension pass" begin
        u = eV(3.0, 2)
        v = MeV(1.0, 2)
        @check_EU_dimension(2, u)
        @check_EU_dimension(2, u, v)
        @test true
    end

    @testset "@check_EU_dimension fail" begin
        u = eV(3.0, 2)
        v = eV(1.0, 1)
        @test_throws DimensionMismatch @check_EU_dimension(1, u)
        @test_throws DimensionMismatch @check_EU_dimension(2, u, v)
    end

    @testset "@check_positive_value pass" begin
        u = eV(3.0, 1)
        v = MeV(2.0, 1)
        @check_positive_value(u)
        @check_positive_value(u, v)
        @test true
    end

    @testset "@check_positive_value fail" begin
        u = eV(-1.0, 1)
        v = eV(1.0, 1)
        @test_throws ArgumentError @check_positive_value(u)
        @test_throws ArgumentError @check_positive_value(v, u)
    end

    @testset "@check_nonnegative_value pass" begin
        u = eV(0.0, 1)
        v = MeV(2.0, 1)
        @check_nonnegative_value(u)
        @check_nonnegative_value(u, v)
        @test true
    end

    @testset "@check_nonnegative_value fail" begin
        u = eV(-1.0, 1)
        v = eV(0.0, 1)
        @test_throws ArgumentError @check_nonnegative_value(u)
        @test_throws ArgumentError @check_nonnegative_value(v, u)
    end
end

@testset "NaturalUnit EUval" begin
    nu = NaturalUnit(GeV)
    u = MeV(1.0, 1)
    val = EUval(nu, u)
    @test val ≈ 1.0e-3
end
