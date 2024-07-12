@with_kw struct ParticleParameter{T<:EnergyUnit} <: MyParameter{T}
    m::EnergyUnit # mass dimension: 1
    α::Real # mass dimension: 0

    unit::Type{T} = T

    @assert m > zero(eV) "m must be positive."
    @assert α > 0 "α must be positive."
end

@with_kw struct ParticleParameterValue{T<:EnergyUnit} <: MyParameterValue{T}
    m::Real
    α::Real

    unit::Type{T} = T
    ndim::Tuple{Int, Int, Int} = (1, 0)

    @assert m > 0 "m must be positive."
    @assert α > 0 "α must be positive."
    @assert ndim == (1, 0) "ndim must be (1, 0, 0)."
end

val(param::ParticleParameter) = val(param.unit, param)
function val(u::Type{<:EnergyUnit}, param::ParticleParameter)
    return ParticleParameterValue{u}(
        val(u, param.m),
        param.α,
        u,
        (1, 0)
    )
end
