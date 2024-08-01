module SimpleScaleFactor

using NaturalUnits
using NonlinearSolve
using OrderedCollections

import SciMLBase: successful_retcode

export set_parameters, scale_factor, Hubble_parameter

const EU = MeV # using MeV as default energy unit.
const NU = NaturalUnit(MeV)

include("parameters.jl")
include("scale_factors.jl")
include("Hubble_parameters.jl")

end # module SimpleScaleFactor
