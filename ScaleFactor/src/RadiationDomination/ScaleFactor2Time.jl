export scale_factor_to_time

function scale_factor_to_time(
    ::Val{:RD}, a::Real;
    T₀=GeV(100), a₀::Real=1, t₀=MeV(0, -1), EU::Type{<:EnergyUnit}=MeV
)
    @check_EU_dimension a 0
    @check_positive_value a

    @check_EU_dimension T₀ 1
    @check_positive_value T₀

    @check_EU_dimension a₀ 0
    @check_positive_value a₀

    @check_EU_dimension t₀ -1

    param = AllParameter()
    param.T₀ = convert(EU, T₀)
    param.a₀ = a₀
    param.EU = EU
    param.NU = NaturalUnit(EU)

    t₀_times_EU = EUval(EU, t₀)
    prob = ODEProblem(__diff_eq_scale_factor_to_time, t₀_times_EU, (a₀, a), param)
    sol = solve(prob)
    @assert successful_retcode(sol)
    return EU(sol(a), -1)
end

function __diff_eq_scale_factor_to_time(t, p, a)
    EU = p.EU
    M_Pl = p.NU.M_Pl
    T = __temperture(a; T₀=p.T₀, a₀=p.a₀, EU=EU)
    H = (π * T^2 / M_Pl) * sqrt(g_star(T) / 90)
    return EUval(EU, 1 / (H * a))
end
