module BlackHoleEvolution

using NaturalUnits

export BH_mass, τ_BH

function BH_mass(M0::EnergyUnit, t::EnergyUnit; t0=MeV(0, -1), g::Real=1)
    @check_EU_dimension M0 1
    @check_EU_dimension t -1
    @check_nonnegative_value M0
    @check_positive_value g

    t < t0 && return M0
    t > t0 && return zero(M0)
    return M0 * cbrt(1 - (t - t0) / τ_BH(M0; g=g))
end

function τ_BH(M0::EnergyUnit; g::Real=1)
    @check_EU_dimension M0 1
    @check_nonnegative_value M0
    @check_positive_value g

    NU = (NaturalUnit ∘ typeof)(M0)
    M_Pl = NU.M_Pl
    return M0^3 / (192 * π^2 * g * M_Pl^4)
end

end # module BlackHoleEvolution
