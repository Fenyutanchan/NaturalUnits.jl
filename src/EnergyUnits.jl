abstract type EnergyUnit end

struct eV <: EnergyUnit
    value # eV(num) means num eV
    dimension::Union{Integer, Rational} # eV(num, dim) means num eV^dim

    eV() = new(1, 1)
    eV(value) = new(value, 1)
    eV(value, dimension) = iszero(dimension) ? value : new(value, dimension)
end

const __head_num_dict = Dict{String, Integer}(
    "k" => 1e3,
    "M" => 1e6,
    "G" => 1e9,
    "T" => 1e12
)

function generation_template_eV(head)
    @assert haskey(__head_num_dict, head)
    return generation_template_eV(head, __head_num_dict[head])
end
function generation_template_eV(head, num)
    heV = Symbol(head, "eV")
    return quote
        export $(heV)

        struct $(heV) <: EnergyUnit
            value # $(head)eV(num) means num $(head)eV
            dimension::Union{Integer, Rational} # $(head)eV(num, dim) means num $(head)eV^dim

            $(heV)() = new(1, 1)
            $(heV)(value) = new(value, 1)
            $(heV)(value, dimension) = iszero(dimension) ? value : new(value, dimension)
        end

        convert(::Type{<:eV}, u::$(heV)) = eV(EUval(u) * $(num)^EUdim(u), EUdim(u))
        convert(::Type{<:$(heV)}, u::eV) = $(heV)(EUval(u) / $(num)^EUdim(u), EUdim(u))
    end
end

for head ∈ keys(__head_num_dict)
    (eval ∘ generation_template_eV)(head)
end

one(u::T) where T<:EnergyUnit = T(1, EUdim(u))
one(::Type{T}) where T<:EnergyUnit = T()

zero(u::T) where T<:EnergyUnit = T(0, EUdim(u))
zero(::Type{T}) where T<:EnergyUnit = T(0) 

convert_EnergyUnit_value_type(T1::Type, u::T2) where T2<:EnergyUnit = T2(convert(T1, EUval(u)), EUdim(u))
convert(::Type{T}, u::T) where T<:EnergyUnit = identity(u)
convert(::Type{<:EnergyUnit}, num::Number) = identity(num)
convert(T::Type{<:EnergyUnit}, u::EnergyUnit) = convert(T, convert(eV, u))

promote_rule(::Type{<:EnergyUnit}, ::Type{<:EnergyUnit}) = eV

__is_same_dimension(u1::EnergyUnit, u2::EnergyUnit) = EUdim(u1) == EUdim(u2)
__diff_dimension_error(u1::EnergyUnit, u2::EnergyUnit, operate::String) = ArgumentError("Cannot $operate energy units with different dimensions: $(EUdim(u1)) and $(EUdim(u2)).")

+(u1::T, u2::T) where T<:EnergyUnit = if __is_same_dimension(u1, u2)
    T(EUval(u1) + EUval(u2), EUdim(u1))
else
    (throw ∘ __diff_dimension_error)(u1, u2, "add")
end
-(u1::T, u2::T) where T<:EnergyUnit = if __is_same_dimension(u1, u2)
    T(EUval(u1) - EUval(u2), EUdim(u1))
else
    (throw ∘ __diff_dimension_error)(u1, u2, "subtract")
end
+(u1::EnergyUnit, u2::EnergyUnit) = +(promote(u1, u2)...)
-(u1::EnergyUnit, u2::EnergyUnit) = -(promote(u1, u2)...)
-(u::T) where T<:EnergyUnit = T(-EUval(u), EUdim(u))

*(u1::T, u2::T) where T<:EnergyUnit = T(EUval(u1) * EUval(u2), EUdim(u1) + EUdim(u2))
/(u1::T, u2::T) where T<:EnergyUnit = T(EUval(u1) / EUval(u2), EUdim(u1) - EUdim(u2))
//(u1::T, u2::T) where T<:EnergyUnit = T(EUval(u1) // EUval(u2), EUdim(u1) - EUdim(u2))

*(u1::EnergyUnit, u2::EnergyUnit) = *(promote(u1, u2)...)
/(u1::EnergyUnit, u2::EnergyUnit) = /(promote(u1, u2)...)
//(u1::EnergyUnit, u2::EnergyUnit) = //(promote(u1, u2)...)

*(num, u::T) where T<:EnergyUnit = T(num * EUval(u), EUdim(u))
/(num, u::T) where T<:EnergyUnit = T(num / EUval(u), -EUdim(u))
//(num, u::T) where T<:EnergyUnit = T(num // EUval(u), -EUdim(u))
*(u::T, num) where T<:EnergyUnit = T(EUval(u) * num, EUdim(u))
/(u::T, num) where T<:EnergyUnit = T(EUval(u) / num, EUdim(u))
//(u::T, num) where T<:EnergyUnit = T(EUval(u) // num, EUdim(u))
^(u::T, num) where T<:EnergyUnit = T(EUval(u)^num, EUdim(u)*num)

==(u1::T, u2::T) where T<:EnergyUnit = EUval(u1) == EUval(u2) && EUdim(u1) == EUdim(u2)
==(u1::EnergyUnit, u2::EnergyUnit) = ==(promote(u1, u2)...)
isless(u1::T, u2::T) where T<:EnergyUnit = if __is_same_dimension(u1, u2)
    isless(EUval(u1), EUval(u2))
else
    (throw ∘ __diff_dimension_error)(u1, u2, "compare")
end
isless(u1::EnergyUnit, u2::EnergyUnit) = isless(promote(u1, u2)...)

sqrt(u::T) where T<:EnergyUnit = T((sqrt ∘ EUval)(u), EUdim(u) // 2)
cbrt(u::T) where T<:EnergyUnit = T((cbrt ∘ EUval)(u), EUdim(u) // 3)

isinf(u::EnergyUnit) = (isinf ∘ EUval)(u) || (isinf ∘ EUdim)(u)
isnan(u::EnergyUnit) = (isnan ∘ EUval)(u) || (isnan ∘ EUdim)(u)

iterate(u::EnergyUnit) = (u, nothing)
iterate(::EnergyUnit, ::Nothing) = nothing
length(u::EnergyUnit) = 1

EUdim(u::EnergyUnit) = u.dimension
EUdim(num::Number) = 0
EUval(u::EnergyUnit) = u.value
EUval(::Type{T}, u::EnergyUnit) where T<:EnergyUnit = convert(T, u).value
EUval(num::Number) = identity(num)
EUval(::Type{T}, num::Number) where T<:EnergyUnit = identity(num)
