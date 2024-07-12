module NaturalUnits

import Base: convert, promote_rule
import Base: +, -, *, /, //, ^, ==, isless
import Base: sqrt, cbrt, one, zero
import Base: getproperty

export EnergyUnit, eV, generation_template_eV
export NaturalUnit
export dim, val

include("NUConstants.jl")
include("SIConstants.jl")

include("MainUnits.jl")
include("NaturalUnit.jl")
include("UnitFunctions.jl")
include("ConstantFunctions.jl")

end # module NaturalUnits
