using AllParameters, NaturalUnits
using CSV, DataFrames, JLD2#=, ProgressMeter=#
using DifferentialEquations#=, =Integrals=#

using SciMLBase: successful_retcode

function main()
    param = AllParameter()
    param.EU = MeV
    param.NU = NaturalUnit(param.EU)
    param.z_eq = 3402
    param.a₀ = 1
    param.a_eq = param.a₀ // (1 + param.z_eq)
    param.t_eq = 51100 * 365 * 24 * 3600 * param.NU.s
    param.t_initial = 1e-11 * param.NU.s
    param.T_eq = 0.8 * one(eV)
    param.t_span = [param.t_eq, param.t_initial]
    
    # progress = Progress(10000)
    # update!(progress, 0)
    # t_times_eV_span = EUval.(Ref(param.EU), param.t_span)
    # @show t_times_eV_span
    # function Friedmann_equation(a, param, t_times_eV)
    #     M_Pl = param.NU.M_Pl
    #     T = T_a(a, param)
    #     @show a, T, a * T
    #     result = π * a * T^2 / (3 * M_Pl) * sqrt(g_star(T) / 10)

    #     ps = (
    #         t_times_eV - first(t_times_eV_span)
    #     ) / (
    #         last(t_times_eV_span) - first(t_times_eV_span)
    #     ) * 10000
    #     update!(progress, floor(Int, ps))

    #     return EUval(param.EU, result)
    # end
    # problem = ODEProblem(
    #     Friedmann_equation,
    #     convert(BigFloat, param.a_eq),
    #     map(BigFloat, t_times_eV_span),
    #     param
    # )
    # sol = solve(problem)
    # @assert successful_retcode(sol)
    # return sol

    log_t_times_eV_span = (log ∘ EUval).(Ref(param.EU), param.t_span)
    # progress = Progress(10000, desc="Solving a-t:")
    # update!(progress, 0)
    function dloga_dlogt(log_a, param, log_t_times_eV)
        a = exp(log_a)
        t = exp(log_t_times_eV) / one(param.EU)
        M_Pl = param.NU.M_Pl
        T = T_a(a, param)
        result = π * t * T^2 / (3 * M_Pl) * sqrt(g_star(T) / 10)

        # ps = (
        #     log_t_times_eV - first(log_t_times_eV_span)
        # ) / (
        #     last(log_t_times_eV_span) - first(log_t_times_eV_span)
        # ) * 10000
        # update!(progress, floor(Int, ps))

        return EUval(param.EU, result)
    end
    prob1 = ODEProblem(
        dloga_dlogt,
        (log ∘ convert)(BigFloat, param.a_eq),
        map(BigFloat, log_t_times_eV_span),
        param
    )
    @info "Solving a-t..."
    sol1 = solve(prob1)
    if successful_retcode(sol1)
        @info "a-t solution is got, saving data..."
        jldopen(joinpath(@__DIR__, "data", "ln_a-ln_t(over_MeV).jld2"), "w") do io
            io["ln(a)"] = sol1.u
            io["ln(t/MeV)"] = sol1.t
        end
        df1 = DataFrame("ln(t/MeV)" => sol1.t, "ln(a)" => sol1.u)
        CSV.write(joinpath(@__DIR__, "data", "ln_a-ln_t(over_MeV).csv"), df1)
    else
        @info "a-t solution is failed to got."
    end
    
    log_a_span = (
        log(param.a_eq),
        last(sol1.u)
    )
    # progress = Progress(10000, desc="Solving t-a:")
    # update!(progress, 0)
    function dlogt_dloga(log_t_times_eV, param, log_a)
        a = exp(log_a)
        t = exp(log_t_times_eV) / one(param.EU)
        M_Pl = param.NU.M_Pl
        T = T_a(a, param)
        result = 3 * M_Pl / (π * t * T^2) * sqrt(10 / g_star(T))

        # ps = (
        #     log_a - first(log_a_span)
        # ) / (
        #     last(log_a_span) - first(log_a_span)
        # ) * 10000
        # update!(progress, floor(Int, ps))

        return EUval(param.EU, result)
    end
    prob2 = ODEProblem(
        dlogt_dloga,
        (log ∘ BigFloat ∘ EUval)(param.EU, param.t_eq),
        map(BigFloat, log_a_span),
        param
    )
    @info "Solving t-a..."
    sol2 = solve(prob2)
    if successful_retcode(sol2)
        @info "t-a solution is got, saving data..."
        jldopen(joinpath(@__DIR__, "data", "ln_t-ln_a(over_MeV).jld2"), "w") do io
            io["ln(t/MeV)"] = sol2.u
            io["ln(a)"] = sol2.t
        end
        df2 = DataFrame("ln(a)" => sol2.t, "ln(t/MeV)" => sol2.u)
        CSV.write(joinpath(@__DIR__, "data", "ln_t-ln_a(over_MeV).csv"), df2)
    else
        @info "t-a solution is failed to got."
    end

    return nothing

    # function integrand(a, param)
    #     M_Pl = param.NU.M_Pl
    #     T = T_a(a, param)
    #     result = 3 * M_Pl / (π * a * T^2) * sqrt(10 / g_star(T))

    #     return EUval(param.EU, result)
    # end
    # pm = Progress(10000, desc="t: ")
    # update!(pm, 0)
    # xpt_step = log(10) / 100
    # a_list = BigFloat[param.a_eq]
    # t_times_eV_list = BigFloat[EUval(param.EU, param.t_eq)]
    # a_terminal = BigFloat(param.a_eq / 10^xpt_step)
    # while true
    #     problem = IntegralProblem(integrand, (first(a_list), a_terminal), param)
    #     sol = solve(problem, QuadGKJL())
    #     @assert successful_retcode(sol)
    #     t_times_eV = sol.u + first(t_times_eV_list)

    #     ps = floor(Int,
    #         EUval(param.EU, sol.u / (param.t_initial - param.t_eq)) * 10000
    #     )
    #     ps < 10001 || break
    #     # update!(pm, ps)

    #     push!(a_list, a_terminal)
    #     push!(t_times_eV_list, t_times_eV)

    #     a_terminal = a_terminal / 10^xpt_step
    #     @show a_terminal, t_times_eV
    # end
end

function g_star(T)
    return 106 + 3 // 4
end
function g_star(EU::Type{<:EnergyUnit}, T)
    return (g_star ∘ EU)(T)
end

function g_star_entropy(T)
    return g_star(T)
end
function g_star_entropy(EU::Type{<:EnergyUnit}, T)
    return (g_star_entropy ∘ EU)(T)
end

function T_a(a, param)
    EU = param.EU
    a_eq = param.a_eq
    T_eq = param.T_eq

    target = g_star_entropy(T_eq) * a_eq^3 * T_eq^3
    target_val = EUval(EU, target)
    # target_dim = EUdim(target)
    f(T_over_EU, EU) = g_star_entropy(EU, T_over_EU) * T_over_EU^3 * a^3 - target_val
    # Tᵢ = cbrt(target / (g_star_entropy(T_eq) * a^3))
    # Tᵢ_over_EU = EUval(EU, Tᵢ)
    T_over_EU_span = EUval.(EU, [T_eq / 2, GeV(1e3)])

    sol = (solve ∘ IntervalNonlinearProblem)(f, T_over_EU_span, EU)
    @assert successful_retcode(sol)
    return EU(sol.u)
end

function a_T(T, param)
    a_eq = param.a_eq
    T_eq = param.T_eq

    target = g_star_entropy(T_eq) * a_eq^3 * T_eq^3
    return cbrt(target / (g_star_entropy(T) * T^3))
end

main()
