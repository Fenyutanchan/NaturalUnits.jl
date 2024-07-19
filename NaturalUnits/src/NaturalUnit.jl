struct NaturalUnit{T<:EnergyUnit}
    unit::Type{T}
end

# in NaturalUnits.jl/src/UnitFunctions.jl 
function __Joule end
function __meter end
function __centimeter end
function __second end
function __kilogram end
function __gram end
function __Kelvin end

# in NaturalUnits.jl/src/ConstantFunctions.jl
function __G_Newton end
function __Planck_MASS end
function __Planck_mass end

__property_function_dict = Dict{Symbol, Function}(
    :J => __Joule,
    :m => __meter,
    :cm => __centimeter,
    :s => __second,
    :kg => __kilogram,
    :g => __gram,
    :K => __Kelvin,

    :G_N => __G_Newton,
    :M_Pl => __Planck_MASS,
    :m_Pl => __Planck_mass
)
function add_property_function(name::Symbol, func::Function)
    __property_function_dict[name] = func
    return nothing
end

function getproperty(u::NaturalUnit, name::Symbol)
    if haskey(__property_function_dict, name)
        return __property_function_dict[name](u)
    else
        return getfield(u, name)
    end
end

EUval(nu::NaturalUnit, x) = EUval(nu.unit, x)
