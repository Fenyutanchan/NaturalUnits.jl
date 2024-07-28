function time_to_scale_factor(
    ::Val{:MR}, t::EnergyUnit;
    T₀=eV(.8), a₀::Real=1/3403, t₀=MeV(0, -1), EU::Type{<:EnergyUnit}=MeV # subscript 0 here means eq.
)
    @check_EU_dimension t -1

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

    t_times_EU_span = map(x -> EUval(EU, x), (t₀, t))
    prob = ODEProblem(__diff_eq_time_to_scale_factor, a₀, t_times_EU_span, param)
    sol = solve(prob)
    @assert successful_retcode(sol)
    return (sol ∘ EUval)(EU, t)
end

function __diff_eq_time_to_scale_factor(a, p, t)
    EU = p.EU
    G_N = p.NU.G_N
    a₀ = p.a₀

    ρₘ = p.ρ₀ₘ * (a₀ / a)^3
    ρᵣ = p.ρ₀ᵣ * (a₀ / a)^4
    ρ = ρₘ + ρᵣ

    H = sqrt(8π * G_N * ρ / 3)

    return EUval(EU, H * a)
end
