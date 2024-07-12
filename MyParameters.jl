module MyParameters

using NaturalUnits
using Parameters

export PrimodialBlackHoleParameter

export change_energy_unit

abstract type MyParameter{T} end

# import NaturalUnits: val

# @with_kw struct DifferentialEquationParameter
# end

# abstract type Parameter{EnergyUnit} end

# struct ParticleParameter
#     m::Union{EnergyUnit, Real}
#     α::Real

#     u::NaturalUnit

#     function ParticleParameter(m, α, u::NaturalUnit=NaturalUnit(MeV))
#         if isa(m, EnergyUnit)
#             @assert dim(m) == 1 "m must be mass unit of 1-dimensional."
#         end
#         return new(m, α, u)
#     end
# end

@with_kw struct PrimodialBlackHoleParameter{T<:EnergyUnit} <: MyParameter{T}
    M₀::Real # mass dimension: 1
    T₀::Real # mass dimension: 1
    β::Real # mass dimension: 0

    type::Type{PrimodialBlackHoleParameter} = PrimodialBlackHoleParameter
    unit::Type{T} = T
    ndim::Tuple{Int, Int, Int} = (1, 1, 0)

    @assert M₀ > 0 "M₀ must be positive."
    @assert T₀ > 0 "T₀ must be positive."
    @assert β > 0 "β must be positive."

    @assert type == PrimodialBlackHoleParameter "Type must be `PrimodialBlackHoleParameter``."
    @assert unit == T "Unit must be the same as $T."
    @assert ndim == (1, 1, 0) "Dimension must be (1, 1, 0)."
end
PrimodialBlackHoleParameter(M₀, T₀, β, unit::Type{T}) where T<:EnergyUnit =
    PrimodialBlackHoleParameter(M₀=M₀, T₀=T₀, β=β, unit=unit)

function change_energy_unit(param::MyParameter, u::Type{<:EnergyUnit})
    u == param.unit && return param
    value_name = fieldnames(param.type)[begin:length(param.ndim)]
    value = Vector{Real}(undef, length(param.ndim))
    for (ii, ndim) ∈ enumerate(param.ndim)
        unitful_value = param.unit(getfield(param, value_name[ii]), ndim)
        value[ii] = (val ∘ convert)(u, unitful_value)
    end
    return param.type(value..., u)
end

end

using NaturalUnits
using .MyParameters

x = PrimodialBlackHoleParameter(1, 2, 3, MeV)
