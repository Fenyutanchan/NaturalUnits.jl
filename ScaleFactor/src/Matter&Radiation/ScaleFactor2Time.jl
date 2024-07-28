function scale_factor_to_time(
    ::Val{:MR}, a::Float64;
    T₀=eV(.8), a₀::Real=1/3403, t₀=MeV(0, -1), EU::Type{<:EnergyUnit}=MeV # subscript 0 here means eq.
)
    @check_EU_dimension a 0
    @check_positive_value a

    @check_EU_dimension T₀ 1
    @check_positive_value T₀

    @check_EU_dimension a₀ 0
    @check_positive_value a₀

    @check_EU_dimension t₀ -1

    param = AllParameter()
    param.EU = EU
    param.NU = NaturalUnit(EU)
    param.a₀ = a₀
    param.ρ₀ₘ = param.ρ₀ᵣ = (π^2 / 30) * g_star(T₀) * T₀^4

    t₀_times_EU = EUval(EU, t₀)
    prob = ODEProblem(__diff_eq_scale_factor_to_time, t₀_times_EU, (a₀, a), param)
    sol = solve(prob)
    @assert successful_retcode(sol)
    return EU(sol(a), -1)
end

function __diff_eq_scale_factor_to_time(t, p, a)
    EU = p.EU
    G_N = p.NU.G_N
    a₀ = p.a₀

    ρₘ = p.ρ₀ₘ * (a₀ / a)^3
    ρᵣ = p.ρ₀ᵣ * (a₀ / a)^4
    ρ = ρₘ + ρᵣ
    H = sqrt(8π * G_N * ρ / 3)

    return EUval(EU, 1 / (H * a))
end
