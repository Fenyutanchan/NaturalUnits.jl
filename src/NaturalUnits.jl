module NaturalUnits

import Base: convert, promote_rule
import Base: +, -, *, /, //, ^, ==, isless
import Base: sqrt, cbrt, one, zero
import Base: isnan, isinf, iszero
import Base: getproperty
import Base: iterate, length # for broadcast

export EnergyUnit, eV, generation_template_eV
export NaturalUnit
export EUdim, EUval
export convert_EnergyUnit_value_type
export add_property_function

include("NUConstants.jl")
include("SIConstants.jl")

include("EnergyUnits.jl")
include("NaturalUnit.jl")
include("UnitFunctions.jl")
include("ConstantFunctions.jl")

end # module NaturalUnits
