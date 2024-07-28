module ScaleFactorRD

using AllParameters
using DifferentialEquations
using EffectiveRelativisticFreedom
using NaturalUnits
using NonlinearSolve

import SciMLBase: successful_retcode

export scale_factor_to_time

function scale_factor_to_time(a::Real; T₀=GeV(100), a₀=1, t₀=MeV(0, -1), EU::Type{<:EnergyUnit}=MeV)
    @check_positive_value a
    @check_positive_value T₀
    @check_positive_value a₀
    @check_nonnegative_value t₀

    param = AllParameter()
    param.T₀ = convert(EU, T₀)
    param.a₀ = a₀
    param.EU = EU
    param.NU = NaturalUnit(EU)

    t₀_times_EU = EUval(EU, t₀)
    prob = ODEProblem(__diff_eq_scale_factor_to_time, t₀_times_EU, (a₀, a), param)
    sol = solve(prob)
    @assert successful_retcode(sol)
    return EU(sol.u[end], -1)
end

function __diff_eq_scale_factor_to_time(t, p, a)
    EU = p.EU
    M_Pl = p.NU.M_Pl
    T = __temperture(a; T₀=p.T₀, a₀=p.a₀, EU=EU)
    H = (π * T^2 / M_Pl) * sqrt(g_star(T) / 90)
    return EUval(EU, 1 / (H * a))
end

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

end # module ScaleFactorRD
