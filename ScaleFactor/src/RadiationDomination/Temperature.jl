function __temperture(a::Real; T₀=GeV(100), a₀=1, EU::Type{<:EnergyUnit}=MeV)
    @check_positive_value a
    @check_positive_value T₀
    @check_positive_value a₀

    param = AllParameter()
    param.EU = EU
    param.target = g_star(T₀) * a₀^3 * T₀^3 / a^3
    
    T_range = (eV(.1), TeV(1.))
    T_over_EU_range = map(x -> EUval(EU, x), T_range)

    prob = IntervalNonlinearProblem(__entropy_eq, T_over_EU_range, param)
    sol = solve(prob)
    @assert successful_retcode(sol)

    return EU(sol.u)
end

function __entropy_eq(T::EnergyUnit, param)
    @check_EU_dimension T 1
    @check_positive_value T

    target = param.target
    @check_EU_dimension target 3
    @check_positive_value target

    return g_star_entropy(T) * T^3 - target
end
function __entropy_eq(T::Real, param)
    EU = param.EU
    return EUval(EU, __entropy_eq(EU(T), param))
end
