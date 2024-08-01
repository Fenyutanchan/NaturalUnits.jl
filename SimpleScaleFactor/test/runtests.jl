using DataFrames
using LaTeXStrings
using NaturalUnits
using Plots
using SimpleScaleFactor
using Test

@testset "parameters" begin
    param = SimpleScaleFactor.set_parameters()
    df = DataFrame(name=String[], value=Any[])
    for (key, value) in param
        push!(df, [key, value])
    end
    @info "Parameters initialized:"
    println(df)

    @testset "time sorting" begin
        t_0 = param["t_0"]
        t_mΛ = param["t_mΛ"]
        t_m = param["t_m"]
        t_eq = param["t_eq"]
        t_r = param["t_r"]
        t_ini = param["t_ini"]
        @test t_0 > t_mΛ > t_m > t_eq > t_r > t_ini
    end

    @info "η-t relation:"
    @testset "η-t function" begin
        df = DataFrame(name=String[], η=EnergyUnit[], η_calc=EnergyUnit[])

        η_eq = param["η_eq"]
        t_eq = param["t_eq"]
        η_eq_calc = SimpleScaleFactor.η_t_fun(t_eq, param)
        push!(df, ["η_eq", η_eq, η_eq_calc])
        @test η_eq / η_eq_calc ≈ 1

        η_m = param["η_m"]
        t_m = param["t_m"]
        η_m_calc = SimpleScaleFactor.η_t_fun(t_m, param)
        push!(df, ["η_m", η_m, η_m_calc])
        @test η_m / η_m_calc ≈ 1

        η_r = param["η_r"]
        t_r = param["t_r"]
        η_r_calc = SimpleScaleFactor.η_t_fun(t_r, param)
        push!(df, ["η_r", η_r, η_r_calc])
        @test η_r / η_r_calc ≈ 1

        for (ii, rr) ∈ (enumerate ∘ sort! ∘ rand)(100)
            η = η_r * (η_m / η_r)^rr
            t = SimpleScaleFactor.t_η_fun(η, param)
            η_calc = SimpleScaleFactor.η_t_fun(t, param)
            push!(df, ["η #$ii", η, η_calc])
            @test η / η_calc ≈ 1
        end

        println(df)
    end
end

@testset "scale factor" begin
    param = SimpleScaleFactor.set_parameters()
    t_ini = param["t_ini"]
    t_r = param["t_r"]
    t_m = param["t_m"]
    t_mΛ = param["t_mΛ"]
    t_0 = param["t_0"]

    df = DataFrame(name=String[], t=EnergyUnit[], scale_factor=Float64[])

    @test begin # before initial time
        push!(df, ["before initial time", t_ini / 2, scale_factor(t_ini / 2)])
        true
    end

    for (ii, rr) ∈ (enumerate ∘ sort! ∘ rand)(100)
        t = t_ini * (t_r / t_ini)^rr
        @test begin
            push!(df, ["radiation domination #$ii", t, scale_factor(t)])
            true
        end
    end

    for (ii, rr) ∈ (enumerate ∘ sort! ∘ rand)(100)
        t = t_r * (t_m / t_r)^rr
        @test begin
            push!(df, ["radiaiton-matter comparable #$ii", t, scale_factor(t)])
            true
        end
    end

    for (ii, rr) ∈ (enumerate ∘ sort! ∘ rand)(100)
        t = t_m * (t_mΛ / t_m)^rr
        @test begin
            push!(df, ["matter domination #$ii", t, scale_factor(t)])
            true
        end
    end

    for (ii, rr) ∈ (enumerate ∘ sort! ∘ rand)(100)
        t = t_mΛ * (t_0 / t_mΛ)^rr
        @test begin
            push!(df, ["matter-dark-energy comparable #$ii", t, scale_factor(t)])
            true
        end
    end

    @test begin # future
        push!(df, ["future", t_0 * 2, scale_factor(t_0 * 2)])
        true
    end

    println(df)

    @test begin
        plot(
            df.t[begin+1:end-1] ./ SimpleScaleFactor.NU.s,
            df.scale_factor[begin+1:end-1],
            label=""
        )
        # plot!(xscale=:log10, yscale=:log10)
        xlabel!(L"$t~[\mathrm{s}]$")
        ylabel!(L"$a$")
        (savefig ∘ joinpath)(@__DIR__, "test_scale_factor.pdf")
        true
    end
end

@testset "Hubble parameter" begin
    param = SimpleScaleFactor.set_parameters()
    t_ini = param["t_ini"]
    t_r = param["t_r"]
    t_m = param["t_m"]
    t_mΛ = param["t_mΛ"]
    t_0 = param["t_0"]

    df = DataFrame(name=String[], t=EnergyUnit[], Hubble_parameter=EnergyUnit[])

    @test begin # before initial time
        push!(df, ["before initial time", t_ini / 2, Hubble_parameter(t_ini / 2)])
        true
    end

    for (ii, rr) ∈ (enumerate ∘ sort! ∘ rand)(100)
        t = t_ini * (t_r / t_ini)^rr
        @test begin
            push!(df, ["radiation domination #$ii", t, Hubble_parameter(t)])
            true
        end
    end

    for (ii, rr) ∈ (enumerate ∘ sort! ∘ rand)(100)
        t = t_r * (t_m / t_r)^rr
        @test begin
            push!(df, ["radiaiton-matter comparable #$ii", t, Hubble_parameter(t)])
            true
        end
    end

    for (ii, rr) ∈ (enumerate ∘ sort! ∘ rand)(100)
        t = t_m * (t_mΛ / t_m)^rr
        @test begin
            push!(df, ["matter domination #$ii", t, Hubble_parameter(t)])
            true
        end
    end

    for (ii, rr) ∈ (enumerate ∘ sort! ∘ rand)(100)
        t = t_mΛ * (t_0 / t_mΛ)^rr
        @test begin
            push!(df, ["matter-dark-energy comparable #$ii", t, Hubble_parameter(t)])
            true
        end
    end

    @test begin # future
        push!(df, ["future", t_0 * 2, Hubble_parameter(t_0 * 2)])
        true
    end

    println(df)

    @test begin
        
        plot(
            df.t[begin+1:end-1] ./ SimpleScaleFactor.NU.s,
            df.Hubble_parameter[begin+1:end-1] ./ eV(2.133e-33),
            label=""
        )
        plot!(xscale=:log10, yscale=:log10)
        xlabel!(L"$t~[\mathrm{s}]$")
        ylabel!(L"$H~[{} \times 100~\mathrm{km}~\mathrm{s}^{-1}~\mathrm{Mpc}^{-1}]$")
        (savefig ∘ joinpath)(@__DIR__, "test_Hubble_parameter.pdf")
        true
    end
end
