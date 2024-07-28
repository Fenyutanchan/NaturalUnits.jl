module ScaleFactorRD

using AllParameters
using DifferentialEquations
using NaturalUnits
using EffectiveRelativisticFreedom

import DifferentialEquations.SciMLBase: successful_retcode

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
    prob = ODEProblem(__diff_eq, t₀_times_EU, (a₀, a), param)
    sol = solve(prob)
    @assert successful_retcode(sol)
    return EU(sol.u[end], -1)
end

function __diff_eq(t, p, a)
    EU = p.EU
    M_Pl = p.NU.M_Pl
    T = __temperture(a; T₀=p.T₀, a₀=p.a₀)
    H = (π * T^2 / M_Pl) * sqrt(g_star(T) / 90)
    return EUval(EU, 1 / (H * a))
end

function __temperture(a::Real; T₀=GeV(100), a₀=1)
    @check_positive_value a
    @check_positive_value T₀
    @check_positive_value a₀

    return (cbrt ∘ g_star_entropy)(T₀) * a₀ * T₀ / a
end

end # module ScaleFactorRD
