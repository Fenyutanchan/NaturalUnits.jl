@with_kw struct PrimodialBlackHoleParameter{T<:EnergyUnit} <: MyParameter{T}
    M₀::EnergyUnit # mass dimension: 1
    T₀::EnergyUnit # mass dimension: 1
    β::Real # mass dimension: 0

    unit::Type{T} = T

    @assert M₀ > zero(eV) "M₀ must be positive."
    @assert T₀ > zero(eV) "T₀ must be positive."
    @assert β > 0 "β must be positive."
end

@with_kw struct PrimodialBlackHoleParameterValue{T<:EnergyUnit} <: MyParameterValue{T}
    M₀::Real
    T₀::Real
    β::Real

    unit::Type{T} = T
    ndim::Tuple{Int, Int, Int} = (1, 1, 0)

    @assert M₀ > 0 "M₀ must be positive."
    @assert T₀ > 0 "T₀ must be positive."
    @assert β > 0 "β must be positive."
    @assert ndim == (1, 1, 0) "ndim must be (1, 1, 0)."
end

val(param::PrimodialBlackHoleParameter) = val(param.unit, param)
function val(u::Type{<:EnergyUnit}, param::PrimodialBlackHoleParameter)
    return PrimodialBlackHoleParameterValue{u}(
        val(u, param.M₀),
        val(u, param.T₀),
        param.β,
        u,
        (1, 1, 0)
    )
end
