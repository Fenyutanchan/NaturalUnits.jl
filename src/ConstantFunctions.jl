function __G_Newton(u::NaturalUnit) # 1 G = ... ()eV^{-2}
    GN = 6.67430e-11 * __meter(u)^3 / (
        __kilogram(u) * __second(u)^2
    )
    return GN
end

function __Planck_MASS(u::NaturalUnit) # 1 M_Planck = ... ()eV
    return 1 / sqrt(8 * π * __G_Newton(u))
end

function __Planck_mass(u::NaturalUnit) # 1 m_Planck = ... ()eV
    return 1 / (sqrt ∘ __G_Newton)(u)
end
