using AllParameters
using JLD2
using LaTeXStrings
using NaturalUnits
using Plots

function main()
    jld2_file_name = joinpath(
        @__DIR__,
        "data",
        "TBD.jld2"
    )

    Eᵢ_over_mᵢ_list = load(jld2_file_name, "Eᵢ_over_mᵢ")
    energy_profile_BE_list = load(jld2_file_name, "energy_profile_BE")
    energy_profile_FD_list = load(jld2_file_name, "energy_profile_FD")
    energy_profile_MB_list = load(jld2_file_name, "energy_profile_MB")
    energy_unit = MeV

    plot(
        Eᵢ_over_mᵢ_list,
        energy_profile_BE_list .* one(energy_unit);
        label="Bose-Einstein Distribution",
        legend=:topright,
        xscale=:log10,
        yscale=:log10,
        title="Energy Profile",
        # minorgrid=true,
        linewith=5,
        color=:blue
    )
    plot!(
        Eᵢ_over_mᵢ_list,
        energy_profile_FD_list .* one(energy_unit);
        label="Fermi-Dirac Distribution",
        linewith=5,
        color=:red
    )
    plot!(
        Eᵢ_over_mᵢ_list,
        energy_profile_MB_list .* one(energy_unit);
        label="Maxwell-Boltzmann Distribution",
        linewith=5,
        color=:green
    )
    xlims!(
        (first ∘ findmin)(Eᵢ_over_mᵢ_list),
        (first ∘ findmax)(Eᵢ_over_mᵢ_list)
    )
    ylims!(1e-65, 1e20)
    xlabel!(L"$E_i / m_i$")
    ylabel!(L"$\frac{\mathrm{d} N_i}{\mathrm{d} E_i}~\left[\mathrm{MeV}^{-1}\right]$")
    savefig("figs/energy_profile.pdf")

    return nothing
end

main()
