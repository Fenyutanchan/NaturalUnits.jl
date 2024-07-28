export time_to_scale_factor

function time_to_scale_factor(t::EnergyUnit; T₀=GeV(100), a₀::Real=1, t₀=MeV(0, -1), EU::Type{<:EnergyUnit}=MeV)
    @check_EU_dimension t -1
    @check_positive_value t

    @check_EU_dimension T₀ 1
    @check_positive_value T₀

    @check_EU_dimension a₀ 0
    @check_positive_value a₀

    @check_EU_dimension t₀ -1
    @check_nonnegative_value t₀

    param = AllParameter()
    param.T₀ = convert(EU, T₀)
    param.a₀ = a₀
    param.EU = EU
    param.NU = NaturalUnit(EU)

    t_times_EU_span = map(x -> EUval(EU, x), (t₀, t))

    prob = ODEProblem(__diff_eq_time_to_scale_factor, a₀, t_times_EU_span, param)
    sol = solve(prob)
    @assert successful_retcode(sol)
    return sol.u[end]
end

function __diff_eq_time_to_scale_factor(a, p, t)
    EU = p.EU
    M_Pl = p.NU.M_Pl
    T = __temperture(a; T₀=p.T₀, a₀=p.a₀, EU=EU)
    H = (π * T^2 / M_Pl) * sqrt(g_star(T) / 90)
    return EUval(EU, H * a)
end
