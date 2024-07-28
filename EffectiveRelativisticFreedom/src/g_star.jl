export g_star

function g_star(T::EnergyUnit; mode=:constant, default::Real=(106+3//4))
    if mode == :constant
        return default
    else
        error("Do not implement this mode yet.")
    end
end

g_star(
    T::Real, EU::Type{<:EnergyUnit};
    mode="constant", default::Real=(106+3//4)
) = g_star(EU(T); mode=mode, default=default)
