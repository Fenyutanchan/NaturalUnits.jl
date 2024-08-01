function __set_a_mΛ!(param_dict)
    a_0 = param_dict["a_0"]
    Ω_m = param_dict["Ω_m"]
    Ω_Λ = param_dict["Ω_Λ"]

    a_mΛ = a_0 * cbrt(Ω_m / Ω_Λ)

    param_dict["a_mΛ"] = a_mΛ
    return a_mΛ
end

function __set_a_eq!(param_dict)
    a_0 = param_dict["a_0"]
    Ω_m = param_dict["Ω_m"]
    Ω_r = param_dict["Ω_r"]

    a_eq = a_0 * Ω_r / Ω_m

    param_dict["a_eq"] = a_eq
    return a_eq
end

function __set_H_0!(param_dict)
    h = param_dict["h"]
    H_0 = convert(MeV, h * eV(2.133e-33))

    param_dict["H_0"] = H_0

    return H_0
end

function __set_t_0!(param_dict)
    H_0 = param_dict["H_0"]
    Ω_m = param_dict["Ω_m"]
    Ω_Λ = param_dict["Ω_Λ"]

    t_0 = 2 * (asinh ∘ sqrt)(Ω_Λ / Ω_m) / (3 * H_0 * sqrt(Ω_Λ))

    param_dict["t_0"] = t_0
    return t_0
end

function __set_t_mΛ!(param_dict)
    EU = param_dict["EU"]
    tol = param_dict["tolerance"]

    t_0 = EUval(EU, param_dict["t_0"])
    function fun(t, tol)
        a_md = scale_factor(Val(:md), t, param_dict)
        a_mΛ = scale_factor(Val(:mΛ), t, param_dict)
        this_tol = abs(a_md - a_mΛ) / (sqrt ∘ abs)(a_md * a_mΛ)
        return this_tol - tol
    end
    prob = NonlinearProblem(fun, t_0, tol)
    sol = solve(prob)
    @assert successful_retcode(sol)
    
    t_mΛ = EU(sol.u, -1)
    param_dict["t_mΛ"] = t_mΛ
    return t_mΛ
end

function __set_η_star!(param_dict)
    a_eq = param_dict["a_eq"]
    t_eq = param_dict["t_eq"]

    η_star = 3 * (2 + sqrt(2)) * t_eq / (2 * a_eq)

    param_dict["η_star"] = η_star
    return η_star
end

function __set_η_eq!(param_dict)
    η_star = param_dict["η_star"]

    η_eq = (sqrt(2) - 1) * η_star

    param_dict["η_eq"] = η_eq
    return η_eq
end

function __set_η_m!(param_dict)
    EU = param_dict["EU"]
    tol = param_dict["tolerance"]

    η_eq = param_dict["η_eq"]
    function fun(η_times_EU, tol)
        η = EU(η_times_EU, -1)
        a_η = a_η_fun(η, param_dict)
        t = t_η_fun(η, param_dict)
        a_md = scale_factor(Val(:md), t, param_dict)

        this_tol = abs(a_η - a_md) / (sqrt ∘ abs)(a_η * a_md)
        return this_tol - tol
    end
    prob = NonlinearProblem(fun, EUval(EU, η_eq), tol)
    sol = solve(prob)
    @assert successful_retcode(sol)

    η_m = EU(sol.u, -1)
    param_dict["η_m"] = η_m
    return η_m
end

function __set_t_m!(param_dict)
    η_m = param_dict["η_m"]

    t_m = t_η_fun(η_m, param_dict)

    param_dict["t_m"] = t_m
    return t_m
end

function __set_η_r!(param_dict)
    EU = param_dict["EU"]
    tol = param_dict["tolerance"]
    H_0 = param_dict["H_0"]
    Ω_r = param_dict["Ω_r"]
    α_r = H_0 * sqrt(Ω_r)

    η_eq = param_dict["η_eq"]
    function fun(η_times_EU, p)
        η = EU(η_times_EU, -1)
        a_rm = a_η_fun(η, param_dict)
        a_rd = α_r * η

        this_tol = abs(a_rm - a_rd) / (sqrt ∘ abs)(a_rm * a_rd)

        return this_tol - tol
    end
    prob = NonlinearProblem(fun, EUval(EU, η_eq), tol)
    sol = solve(prob)
    @assert successful_retcode(sol)

    η_r = EU(sol.u, -1)
    param_dict["η_r"] = η_r
    return η_r
end

function __set_t_r!(param_dict)
    η_r = param_dict["η_r"]

    t_r = t_η_fun(η_r, param_dict)

    param_dict["t_r"] = t_r
    return t_r
end

function __set_t_ini!(param_dict)
    NU = param_dict["NU"]
    g_star = param_dict["g_star"]
    T_ini = param_dict["T_ini"]
    M_Pl = NU.M_Pl

    H_sqr = (π^2 / 90) * g_star * T_ini^4 / M_Pl^2

    t_ini = 1 / sqrt(2 * H_sqr)
    param_dict["t_ini"] = t_ini
    return t_ini
end
