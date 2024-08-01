Hubble_parameter(t) = Hubble_parameter(t, set_parameters())

function Hubble_parameter(t::EnergyUnit, param_dict::Dict)
    t_ini = param_dict["t_ini"]
    t_r = param_dict["t_r"]
    t_m = param_dict["t_m"]
    t_mΛ = param_dict["t_mΛ"]

    @check_nonnegative_value t

    if t < t_r
        t < t_ini && @warn "Time is before the initial time!"
        return Hubble_parameter(Val(:rd), t, param_dict)
    elseif t < t_m
        return Hubble_parameter(Val(:rm), t, param_dict)
    elseif t < t_mΛ
        return Hubble_parameter(Val(:md), t, param_dict)
    else
        return Hubble_parameter(Val(:mΛ), t, param_dict)
    end
end
Hubble_parameter(t::Real, param_dict::Dict) = Hubble_parameter(EU(t, -1), param_dict)

Hubble_parameter(::Val{:rd}, t::EnergyUnit, param_dict::Dict) = 1 / (2 * t)

function Hubble_parameter(::Val{:rm}, t::EnergyUnit, param_dict::Dict)
    a_eq = param_dict["a_eq"]
    η_star = param_dict["η_star"]

    η = η_t_fun(t, param_dict)

    return 2 * η_star^2 * (η + η_star) / (
        a_eq * η^2 * (η + 2 * η_star)^2
    )
end

Hubble_parameter(::Val{:md}, t::EnergyUnit, param_dict::Dict) = 2 / (3 * t)

function Hubble_parameter(::Val{:mΛ}, t::EnergyUnit, param_dict::Dict)
    H_0 = param_dict["H_0"]
    Ω_Λ = param_dict["Ω_Λ"]

    α = H_0 * sqrt(Ω_Λ)

    return α * coth((3/2) * α * t)
end

Hubble_parameter(T::Val, t::Real, param_dict::Dict) =
    Hubble_parameter(T, EU(t, -1), param_dict)
    