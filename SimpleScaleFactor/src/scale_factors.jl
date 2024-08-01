scale_factor(t) = scale_factor(t, set_parameters())

function scale_factor(t::EnergyUnit, param_dict::Dict)
    t_ini = param_dict["t_ini"]
    t_r = param_dict["t_r"]
    t_m = param_dict["t_m"]
    t_mΛ = param_dict["t_mΛ"]

    @check_nonnegative_value t

    if t < t_r
        t < t_ini && @warn "Time is before the initial time!"
        return scale_factor(Val(:rd), t, param_dict)
    elseif t < t_m
        return scale_factor(Val(:rm), t, param_dict)
    elseif t < t_mΛ
        return scale_factor(Val(:md), t, param_dict)
    else
        return scale_factor(Val(:mΛ), t, param_dict)
    end
end

function scale_factor(t::Real, param_dict::Dict)
    EU = param_dict["EU"]
    return scale_factor(EU(t, -1), param_dict)
end

function scale_factor(::Val{:rd}, t::EnergyUnit, param_dict::Dict)
    H_0 = param_dict["H_0"]
    Ω_r = param_dict["Ω_r"]

    α = H_0 * sqrt(Ω_r)

    return sqrt(α * t)
end

function scale_factor(::Val{:rm}, t::EnergyUnit, param_dict::Dict)
    η = η_t_fun(t, param_dict)
    return a_η_fun(η, param_dict)
end

function scale_factor(::Val{:md}, t::EnergyUnit, param_dict::Dict)
    H_0 = param_dict["H_0"]
    Ω_m = param_dict["Ω_m"]
    
    α = H_0 * sqrt(Ω_m)

    return cbrt((3/2) * α * t)^2
end

function scale_factor(::Val{:mΛ}, t::EnergyUnit, param_dict::Dict)
    H_0 = param_dict["H_0"]
    Ω_m = param_dict["Ω_m"]
    Ω_Λ = param_dict["Ω_Λ"]

    α = H_0 * sqrt(Ω_Λ)
    return cbrt((Ω_m / Ω_Λ) * sinh((3/2) * α * t)^2)
end

function scale_factor(T::Val, t::Real, param_dict::Dict)
    EU = param_dict["EU"]
    return scale_factor(T, EU(t, -1), param_dict)
end

function a_η_fun(η::EnergyUnit, param_dict)
    a_eq = param_dict["a_eq"]
    η_star = param_dict["η_star"]

    return a_eq * ((η / η_star)^2 + 2 * η / η_star)
end

function a_η_fun(η::Real, param_dict)
    EU = param_dict["EU"]
    return a_η_fun(EU(η, -1), param_dict)
end

function t_η_fun(η::EnergyUnit, param_dict)
    a_eq = param_dict["a_eq"]
    η_star = param_dict["η_star"]

    return a_eq * η^2 * (η + 3 * η_star) / (3 * η_star^2)
end

function t_η_fun(η::Real, param_dict)
    EU = param_dict["EU"]
    return EUval(EU, t_η_fun(EU(η, -1), param_dict))
end

function η_t_fun(t::Real, param_dict)
    EU = param_dict["EU"]
    a_eq = param_dict["a_eq"]
    η_star = param_dict["η_star"]
    η_m = param_dict["η_m"]
    η_r = param_dict["η_r"]

    a = EUval(EU, a_eq / (3 * η_star^2))
    b = EUval(EU, a_eq / η_star)
    c = 0
    d = EUval(EU, -t)

    η_span = map(η -> EUval(EU, η), (η_r/2, 2*η_m))

    fun = (x, p) -> a * x^3 + b * x^2 + c * x + d
    prob = IntervalNonlinearProblem(fun, η_span)
    sol = solve(prob)
    # @show sol.retcode
    # @assert successful_retcode(sol)
    return sol.u

    # init_epsilon = 1e-20
    # while true
    #     root_list = roots([d, c, b, a]; polish=true, epsilon=init_epsilon)
    #     # @show root_list
    #     filter!(isreal, root_list)
    #     if isempty(root_list)
    #         init_epsilon /= 10
    #         continue
    #     end

    #     real_root_list = map(real, root_list)
    #     sort!(real_root_list)
    #     init_epsilon /= 10
    #     return EU(real_root_list[end], -1)
    # end
end

function η_t_fun(t::EnergyUnit, param_dict)
    EU = param_dict["EU"]
    return EU(η_t_fun(EUval(EU, t), param_dict), -1)
end
