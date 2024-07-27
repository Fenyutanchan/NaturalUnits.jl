module BlackHoleEvolution

using NaturalUnits

export BH_mass

function BH_mass(M0::EnergyUnit, t::EnergyUnit; t0=MeV(0, -1), g::Real=1)
    @assert EUdim(M0) == 1 "M0 should have mass dimension of 1, but got $(EUdim(M0))."
    @assert EUdim(t) == -1 "t should have mass dimension of -1, but got $(EUdim(t))."
    @assert M0 > zero(M0) "M0 should be positive, but got $(M0)."
    @assert g > 0 "g should be positive, but got $(g)."

    t < t0 && return M0
    t > t0 && return zero(M0)
    return M0 * cbrt(1 - (t - t0) / τ_BH(M0; g=g))
end

function τ_BH(M0::EnergyUnit; g::Real=1)
    @assert EUdim(M0) == 1 "M0 should have mass dimension of 1, but got $(EUdim(M0))."
    @assert M0 > zero(M0) "M0 should be positive, but got $(M0)."
    @assert g > 0 "g should be positive, but got $(g)."

    NU = (NaturalUnit ∘ typeof)(M0)
    M_Pl = NU.M_Pl
    return M0^3 / (192 * π^2 * g * M_Pl^4)
end

end # module BlackHoleEvolution
