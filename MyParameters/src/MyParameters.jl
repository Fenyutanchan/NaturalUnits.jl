module MyParameters

using NaturalUnits
using Parameters

import NaturalUnits: val

export PrimodialBlackHoleParameter

include("MyParameter.jl")
include("PrimodialBlackHoleParameter.jl")

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

end
