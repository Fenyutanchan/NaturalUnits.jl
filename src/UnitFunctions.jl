function __Joule(u::NaturalUnit) # 1 J = ... ()eV
    J = convert(u.unit, one(eV)) / SIConstants.e
    return convert(u.unit, J)
end

function __meter(u::NaturalUnit) # 1 m = ... ()eV^{-1}
    m = (
        NUConstants.ħ * NUConstants.c / __Joule(u)
    ) / (
        SIConstants.ħ * SIConstants.c
    )
    return m
end
__centimeter(u::NaturalUnit) = __meter(u) * 1e-2 # 1 cm = ... ()eV^{-1}

function __second(u::NaturalUnit) # 1 s = ... ()eV^{-1}
    return SIConstants.c * __meter(u)
end

function __kilogram(u::NaturalUnit) # 1 kg = ... ()eV
    kg = __Joule(u)/ (
        __meter(u)^2 / __second(u)^2
    )
    return kg
end
__gram(u::NaturalUnit) = __kilogram(u) * 1e-3 # 1 g = ... ()eV

function __Kelvin(u::NaturalUnit) # 1 K = ... ()eV
    K = SIConstants.k_B * __Joule(u)
    return K
end
