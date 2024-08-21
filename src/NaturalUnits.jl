# Copyright (c) 2024 Quan-feng WU <wuquanfeng@ihep.ac.cn>
# 
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

module NaturalUnits

import Base: convert, promote_rule
import Base: +, -, *, /, //, ^, ==, inv, isless
import Base: sqrt, cbrt, one, zero
import Base: isnan, isinf, iszero
import Base: getproperty
import Base: iterate, length # for broadcast

export EnergyUnit, eV, generation_template_eV
export NaturalUnit
export EUdim, EUval
export convert_EnergyUnit_value_type
export add_property_function

export @check_EU_dimension, @check_positive_value, @check_nonnegative_value

include("NUConstants.jl")
include("SIConstants.jl")

include("EnergyUnits.jl")
include("NaturalUnit.jl")
include("UnitFunctions.jl")
include("ConstantFunctions.jl")

include("Utilities.jl")

end # module NaturalUnits
