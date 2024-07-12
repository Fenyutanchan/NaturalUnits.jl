module MyParameters

using NaturalUnits
using Parameters

import NaturalUnits: val

export PrimodialBlackHoleParameter

include("MyParameter.jl")
include("ParticleParameter.jl")
include("PrimodialBlackHoleParameter.jl")

# @with_kw struct DifferentialEquationParameter
# end


end
