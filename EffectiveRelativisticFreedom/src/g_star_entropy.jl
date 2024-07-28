export g_star_entropy

function g_star_entropy(T::EnergyUnit; mode=:constant, default::Real=(106+3//4))
    if mode == :constant
        return default
    else
        error("Do not implement this mode yet.")
    end
end

g_star_entropy(
    T::Real, EU::Type{<:EnergyUnit};
    mode="constant", default::Real=(106+3//4)
) = g_star(EU(T); mode=mode; default=default)
