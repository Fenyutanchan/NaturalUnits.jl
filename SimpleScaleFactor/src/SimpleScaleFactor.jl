module SimpleScaleFactor

using NaturalUnits
using NonlinearSolve
using OrderedCollections
using RealPolynomialRoots

import SciMLBase: successful_retcode

export set_parameters, scale_factor

const EU = MeV # using MeV as default energy unit.
const NU = NaturalUnit(MeV)

include("parameters.jl")
include("scale_factors.jl")

end # module SimpleScaleFactor
