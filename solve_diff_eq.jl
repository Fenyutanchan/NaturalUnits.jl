using AllParameters
using DifferentialEquations
using Integrals
using JLD2
using NaturalUnits
using Plots
using ProgressMeter

function main()
    unit = NaturalUnit(MeV)
    param = AllParameter()
    param.unit = unit

    param.g_star = 106.75

    param.M₀_BH = 10 * unit.g
    param.α_BH = .2
    param.β_BH = 1e-6

    param.g_X = 2
    param.m_X = GeV(1e12)
    param.α_X = .5^2 / (16 * π)

    param.a₀ = 1
    param.a_init = 1
    param.na³_init = 0
    param.a_finl = 1.5e6

    param.T₀ = (1 / 2) * (
        5 / (param.g_star * π^3)
    )^(1/4) * sqrt(3 * param.α_BH * unit.m_P^3 / param.M₀_BH)
    param.C_aT = param.T₀ * param.a₀
    param.g_BH = 7e-5 * param.g_star
    param.ρ₀ = (π^2 / 30) * param.g_star * param.T₀^4

    progress = Progress(10000, desc="a: ")
    update!(progress, 0)

    function Boltzmann_function(na³, p, a)
        na³_times_eV = param.unit.unit(na³, 3)
        T = param.C_aT / a
        M_P = param.unit.M_P
        g_star = param.g_star
        inv_LHS_factor = a^2 * sqrt(
            90 * M_P^2 / (π^2 * g_star * T^4)
        )

        ps = (a - param.a_init) / (param.a_finl - param.a_init) * 10000
        update!(progress, ceil(Int, ps))

        RHS = collision_term_production(a, param) - collision_term_depletion(na³_times_eV, a, param)
        return val(param.unit.unit, inv_LHS_factor * RHS)
    end
    problem = ODEProblem(Boltzmann_function, param.na³_init, (param.a_init, param.a_finl))

    solution = solve(problem)

    a_list = solution.t
    na³_list = [param.unit.unit(na³, 3) for na³ in solution.u]

    current_jld2_file_name = joinpath(
        @__DIR__,
        "data",
        "current.jld2"
    )
    jld2_file_name = joinpath(
        @__DIR__,
        "data",
        get_jld2_file_name(param)
    )
    jldopen(jld2_file_name, "w") do file
        file["a"] = a_list
        file["na³"] = na³_list
        file["param"] = param.content
    end
    rm(current_jld2_file_name; force=true, recursive=true)
    symlink(jld2_file_name, current_jld2_file_name)

    return jld2_file_name
end

function get_jld2_file_name(param)
    M₀_BH = param.M₀_BH / param.unit.g
    β_BH = param.β_BH
    m_X = param.m_X
    α_X = param.α_X

    return "M₀_BH_gram($(M₀_BH), 1)_β_BH_$(β_BH)_m_X_$(m_X)_α_X_$(α_X).jld2"
end

function particle_production_rate_per_time_and_energy(energy, T_BH, param) # mass dimension: 0
    iszero(T_BH) && return zero(param.m_X)
    Γ_X = greybody_factor(energy, T_BH, param)
    g_X = param.g_X

    return (g_X / (2 * π)) * Γ_X * exp(-energy / T_BH)
end

function greybody_factor(energy, T_BH, param) # mass dimension: 0
    iszero(T_BH) && return 0
    p_X_sqr = energy^2 - param.m_X^2
    p_X_sqr ≤ zero(p_X_sqr) && return 0
    return 27 * p_X_sqr / (64 * π^2 * T_BH^2)
end

function particle_production_rate_from_BH(T_BH, param)
    iszero(T_BH) && return zero(param.m_X)

    g_X = param.g_X
    m_X = param.m_X
    prefactor = 27 * g_X / (64 * π^3)

    return prefactor * (m_X + T_BH) * exp(-m_X / T_BH)
end

temperature_BH(M_BH, unit) = M_BH ≤ zero(M_BH) ? zero(M_BH) :
    1 / (8 * π * unit.G_N * M_BH)

function particle_full_number_energy_profile(energy, param) # mass dimension: -1
    τ_BH = val(param.unit.unit, life_BH(param))

    function integrand(t, p)
        t_over_eV = param.unit.unit(t, -1)
        m_BH = mass_BH_evolution_time(t_over_eV, param)
        T_BH = temperature_BH(m_BH, param.unit)
        result = particle_production_rate_per_time_and_energy(energy, T_BH, param)
        # @assert (iszero ∘ dim)(result) "The result must have mass dimension 0."

        return val(param.unit.unit, result)
    end

    problem = IntegralProblem(integrand, (zero(τ_BH), τ_BH))
    integral = solve(problem, QuadGKJL())

    @assert integral.retcode == ReturnCode.Success "Integral calculation failed."

    return param.unit.unit(integral.u, -1)
end

function life_BH(param) # mass dimension: -1
    M₀_BH = param.M₀_BH
    g_BH = param.g_BH
    m_P = param.unit.m_P
    return M₀_BH^3 / (3 * g_BH * m_P^4)
end

function mass_BH_evolution_time(t, param)
    M₀_BH = param.M₀_BH
    τ_BH = life_BH(param)

    t ≤ zero(t) && return M₀_BH
    t ≥ τ_BH && return zero(M₀_BH)
    return M₀_BH * cbrt(1 - t / τ_BH)
end

function mass_BH_evolution(a, param)
    M_P = param.unit.M_P
    g_star = param.g_star
    C_aT = param.C_aT
    a₀ = param.a₀

    prefactor = sqrt(10 / g_star) * 3 * M_P / (π * C_aT^2)
    t = prefactor * (a^2 - a₀^2) / 2

    return mass_BH_evolution_time(t, param)
end

function number_density_BH(a, param)
    M₀_BH = param.M₀_BH
    β_BH = param.β_BH
    ρ₀ = param.ρ₀
    a₀ = param.a₀

    prefactor = β_BH * ρ₀ / M₀_BH
    return prefactor * (a₀ / a)^3
end

function collision_term_production(a, param)
    m_BH = mass_BH_evolution(a, param)
    T_BH = temperature_BH(m_BH, param.unit)
    Γ = particle_production_rate_from_BH(T_BH, param)
    n_BH = number_density_BH(a, param)

    return Γ * n_BH
end

function collision_term_depletion(na³, a, param; mode="quick")
    inv_γ = if mode == "slow"
        error("Not implemented.")
    elseif mode == "quick"
        inv_γ_quick(a, param)
    elseif mode == "quick quick"
        error("Not implemented.")
    else
        error("Invalid mode: $mode.")
    end

    return (na³ / a^3) * particle_decay_rate(inv_γ, param)
end

particle_decay_rate(inv_γ, param) = inv_γ * param.α_X * param.m_X

function inv_γ_quick(a, param)
    global inv_γ_cache
    
    m_X = param.m_X
    m_BH = mass_BH_evolution(a, param)
    T_BH = temperature_BH(m_BH, param.unit)
    ∞ = 709.782712893384 * T_BH

    function num_integrand(E, p)
        E_times_eV = param.unit.unit(E, 1)
        den = particle_production_rate_per_time_and_energy(E_times_eV, T_BH, param)
        num = den * param.m_X / E_times_eV
        return val(param.unit.unit, num)
    end
    function den_integrand(E, p)
        E_times_eV = param.unit.unit(E, 1)
        den = particle_production_rate_per_time_and_energy(E_times_eV, T_BH, param)
        return val(param.unit.unit, den)
    end

    m_X_over_eV = val(param.unit.unit, m_X)
    ∞_over_eV = val(param.unit.unit, ∞)
    integration_region = (m_X_over_eV, ∞_over_eV)

    num_problem = IntegralProblem(num_integrand, integration_region)
    den_problem = IntegralProblem(den_integrand, integration_region)

    num_integral = solve(num_problem, QuadGKJL())
    den_integral = solve(den_problem, QuadGKJL())

    @assert num_integral.retcode == ReturnCode.Success "Integral calculation failed."
    @assert den_integral.retcode == ReturnCode.Success "Integral calculation failed."

    inv_γ = num_integral.u / den_integral.u
    isnan(inv_γ) && return inv_γ_cache
    inv_γ_cache = inv_γ
    return inv_γ
end

main()
