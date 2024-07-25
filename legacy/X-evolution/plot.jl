using JLD2
using LaTeXStrings
using NaturalUnits
using Plots
using SpecialFunctions

function main()
    jld2_file_name = joinpath(
        @__DIR__,
        "data",
        "current.jld2"
    )
    a_list = load(jld2_file_name, "a")
    na³_list = load(jld2_file_name, "na³")
    param_dict = load(jld2_file_name, "param")
    mode = param_dict[:mode]

    ζ₃ = zeta(3)
    T_list = param_dict[:C_aT] ./ a_list
    m_over_T_list = param_dict[:m_X] ./ T_list
    nX_list = na³_list ./ a_list.^3
    nγ_list = 2 * ζ₃ * T_list.^3 / π^2
    nX_over_nγ_list = nX_list ./ nγ_list
    plateau_list = [1.9e-13 for _ in m_over_T_list]

    plot(
        m_over_T_list, nX_over_nγ_list;
        linewith=3,
        label=L"$n_X / n_γ$", legend=:topleft,
        xscale=:log10, yscale=:log10,
        title=L"$(n_X / n_γ)$ - $(m_X / T)$ in mode of " * mode,
        minorgrid=true
    )
    plot!(
        m_over_T_list, plateau_list;
        linewith=3, label=L"$1.9 \times 10^{-13}$",
        color=:red, linestyle=:dash
    )
    xlabel!(L"$m_X / T$")
    ylabel!(L"$n_X / n_γ$")
    xlims!(1e-2, 1.2e3)
    ylims!(1e-14, 1e-7)
    savefig("figs/nX_over_nγ-m_over_T.pdf")
end

main()
