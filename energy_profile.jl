using AllParameters
using NaturalUnits
using Nemo
using JLD2

import Base: convert

RR = RealField()
set_precision!(RR, 1000)

convert(::Type{RealFieldElem}, num::Real) = RR(num)

function main()
    param = AllParameter()

    param.unit_to_storage = MeV
    param.unit = NaturalUnit(param.unit_to_storage)
    
    param.mᵢ = eV(RR(1))
    param.gᵢ = RR(1)
    param.gStar = RR(106.75)
    param.gStarPBH = 7e-5 * param.gStar
    param.mPBH₀ = RR(1) * param.unit.g

    param.α = RR(.2)
    param.β = 1.1e-6 * (param.α / .2)^(-1/2) * (param.mPBH₀ / (1e4 * param.unit.g))^(-17/24)
    param.T₀ = 2 * (10 / param.gStar)^(1/4) * sqrt(3 * param.α * param.unit.M_Pl^3 / param.mPBH₀)
    param.ρR₀ = (const_pi(RR)^2 / 30) * param.gStar * param.T₀^4
    param.ρPBH₀ = param.β * param.ρR₀
    param.nPBH₀ = param.ρPBH₀ / param.mPBH₀
    param.H₀ = sqrt(param.ρR₀ / (3 * param.unit.M_Pl^2))

    @assert (!is_nonzero ∘ EUval)(param.H₀^2 - (const_pi(RR)^2 / 90) * param.gStar * param.T₀^4 / param.unit.M_Pl^2)
    @assert (!is_nonzero ∘ EUval)(param.H₀^2 - param.ρR₀ / (3 * param.unit.M_Pl^2))

    param.exponent_list = 0:.01:40
    jldopen(joinpath(@__DIR__, "data", generate_jld_file_name(param)), "w") do io
        Eᵢ_list = Vector{param.unit_to_storage}(undef, length(param.exponent_list))
        energy_profile_BE_list = deepcopy(Eᵢ_list)
        energy_profile_FD_list = deepcopy(Eᵢ_list)
        energy_profile_MB_list = deepcopy(Eᵢ_list)

        Threads.@threads for ii ∈ eachindex(param.exponent_list)
            Eᵢ = param.mᵢ * 10^param.exponent_list[ii]
            EP_BE = energy_profile_BE(Eᵢ, param)
            EP_FD = energy_profile_FD(Eᵢ, param)
            EP_MB = energy_profile_MB(Eᵢ, param)

            Eᵢ_list[ii] = convert(param.unit_to_storage, Eᵢ)
            energy_profile_BE_list[ii] = convert(param.unit_to_storage, EP_BE)
            energy_profile_FD_list[ii] = convert(param.unit_to_storage, EP_FD)
            energy_profile_MB_list[ii] = convert(param.unit_to_storage, EP_MB)
        end

        io["Eᵢ_over_mᵢ"] = 10 .^ param.exponent_list
        io["energy_profile_BE"] = energy_profile_BE_list
        io["energy_profile_FD"] = energy_profile_FD_list
        io["energy_profile_MB"] = energy_profile_MB_list
    end

    return nothing
end

function generate_jld_file_name(param)
    
    return "TBD.jld2"
end

function energy_profile_BE(Eᵢ, param)
    mᵢ = param.mᵢ
    Eᵢ < mᵢ && return zero(Eᵢ^-1)

    gᵢ = param.gᵢ
    gStarPBH = param.gStarPBH
    mPBH₀ = param.mPBH₀
    M_Pl = param.unit.M_Pl

    exp_expr = exp(-Eᵢ * mPBH₀ / M_Pl^2)

    prefactor = 27 * gᵢ * (Eᵢ^2 - mᵢ^2) / (128 * const_pi(RR)^3)
    integral = (
        1 / (
            64 * Eᵢ^5 * gStarPBH * M_Pl^6 * const_pi(RR)^2
        )
    ) * (
        Eᵢ^4 * mPBH₀^4 * log(1 - exp_expr)
        - 4 * Eᵢ^3 * mPBH₀^3 * M_Pl^2 * polylog(2, exp_expr)
        - 12 * Eᵢ^2 * mPBH₀^2 * M_Pl^4 * polylog(3, exp_expr)
        - 24 * Eᵢ * mPBH₀ * M_Pl^6 * polylog(4, exp_expr)
        - 24 * M_Pl^8 * polylog(5, exp_expr)
        + 24 * M_Pl^8 * zeta(5, RR)
    )

    return prefactor * integral
end

function energy_profile_FD(Eᵢ, param)
    mᵢ = param.mᵢ
    Eᵢ < mᵢ && return zero(Eᵢ^-1)

    gᵢ = param.gᵢ
    gStarPBH = param.gStarPBH
    mPBH₀ = param.mPBH₀
    M_Pl = param.unit.M_Pl

    exp_expr = exp(-Eᵢ * mPBH₀ / M_Pl^2)

    prefactor = 27 * gᵢ * (Eᵢ^2 - mᵢ^2) / (128 * const_pi(RR)^3)
    integral = (
        1 / (
            128 * Eᵢ^5 * gStarPBH * M_Pl^6 * const_pi(RR)^2
        )
    ) * (
        -2 * Eᵢ^4 * mPBH₀^4 * log(1 + exp_expr)
        + 8 * Eᵢ^3 * mPBH₀^3 * M_Pl^2 * polylog(2, -exp_expr)
        + 24 * Eᵢ^2 * mPBH₀^2 * M_Pl^4 * polylog(3, -exp_expr)
        + 48 * Eᵢ * mPBH₀ * M_Pl^6 * polylog(4, -exp_expr)
        + 48 * M_Pl^8 * polylog(5, -exp_expr)
        + 45 * M_Pl^8 * zeta(5, RR)
    )

    return prefactor * integral
end

function energy_profile_MB(Eᵢ, param)
    mᵢ = param.mᵢ
    Eᵢ < mᵢ && return zero(Eᵢ^-1)

    gᵢ = param.gᵢ
    gStarPBH = param.gStarPBH
    mPBH₀ = param.mPBH₀
    M_Pl = param.unit.M_Pl

    prefactor = 27 * gᵢ * (Eᵢ^2 - mᵢ^2) / (128 * const_pi(RR)^3)
    integral = -M_Pl^2 * (
        -24 + gamma(RR(5), Eᵢ * mPBH₀ / M_Pl^2)
    ) / (
        64 * Eᵢ^5 * gStarPBH * const_pi(RR)^2
    )

    return prefactor * integral
end

main()
