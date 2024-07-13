module NaturalUnits

import Base: convert, promote_rule
import Base: +, -, *, /, //, ^, ==, isless
import Base: sqrt, cbrt, one, zero
import Base: isnan, isinf
import Base: getproperty
import Base: iterate, length # for broadcast

export EnergyUnit, eV, generation_template_eV
export NaturalUnit
export dim, val
export add_property_function

include("NUConstants.jl")
include("SIConstants.jl")

include("MainUnits.jl")
include("NaturalUnit.jl")
include("UnitFunctions.jl")
include("ConstantFunctions.jl")

end # module NaturalUnits
