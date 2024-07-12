abstract type EnergyUnit end

struct eV <: EnergyUnit
    value::Real # eV(num) means num eV
    dimension::Union{Integer, Rational} # eV(num, dim) means num eV^dim

    eV() = new(1, 1)
    eV(value) = new(value, 1)
    eV(value, dimension) = iszero(dimension) ? value : new(value, dimension)
end

const __head_num_dict = Dict{String, Int}(
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
        export $heV

        struct $heV <: EnergyUnit
            value::Real # $(head)eV(num) means num $(head)eV
            dimension::Union{Integer, Rational} # $(head)eV(num, dim) means num $(head)eV^dim

            $heV() = new(1, 1)
            $(heV)(value) = new(value, 1)
            $(heV)(value, dimension) = iszero(dimension) ? value : new(value, dimension)
        end

        convert(::Type{eV}, u::$heV) = eV(val(u) * $(num)^dim(u), dim(u))
        convert(::Type{$heV}, u::eV) = $heV(val(u) / $(num)^dim(u), dim(u))
    end
end

for head ∈ keys(__head_num_dict)
    (eval ∘ generation_template_eV)(head)
end

one(::T) where T<:EnergyUnit = T()
one(::Type{T}) where T<:EnergyUnit = T()

zero(::T) where T<:EnergyUnit = T(0)
zero(::Type{T}) where T<:EnergyUnit = T(0) 

convert(::Type{T}, u::T) where T<:EnergyUnit = identity(u)
convert(T::Type{<:EnergyUnit}, u::EnergyUnit) = convert(T, convert(eV, u))

promote_rule(::Type{<:EnergyUnit}, ::Type{<:EnergyUnit}) = eV

__is_same_dimension(u1::EnergyUnit, u2::EnergyUnit) = dim(u1) == dim(u2)
__diff_dimension_error(u1::EnergyUnit, u2::EnergyUnit, operate::String) = ArgumentError("Cannot $operate energy units with different dimensions: $(dim(u1)) and $(dim(u2)).")

+(u1::T, u2::T) where T<:EnergyUnit = if __is_same_dimension(u1, u2)
    T(val(u1) + val(u2), dim(u1))
else
    (throw ∘ __diff_dimension_error)(u1, u2, "add")
end
-(u1::T, u2::T) where T<:EnergyUnit = if __is_same_dimension(u1, u2)
    T(val(u1) - val(u2), dim(u1))
else
    (throw ∘ __diff_dimension_error)(u1, u2, "subtract")
end
+(u1::EnergyUnit, u2::EnergyUnit) = +(promote(u1, u2)...)
-(u1::EnergyUnit, u2::EnergyUnit) = -(promote(u1, u2)...)

*(u1::T, u2::T) where T<:EnergyUnit = T(val(u1) * val(u2), dim(u1) + dim(u2))
/(u1::T, u2::T) where T<:EnergyUnit = T(val(u1) / val(u2), dim(u1) - dim(u2))
//(u1::T, u2::T) where T<:EnergyUnit = T(val(u1) // val(u2), dim(u1) - dim(u2))

*(u1::EnergyUnit, u2::EnergyUnit) = *(promote(u1, u2)...)
/(u1::EnergyUnit, u2::EnergyUnit) = /(promote(u1, u2)...)
//(u1::EnergyUnit, u2::EnergyUnit) = //(promote(u1, u2)...)

*(num, u::T) where T<:EnergyUnit = T(num * val(u), dim(u))
/(num, u::T) where T<:EnergyUnit = T(num / val(u), -dim(u))
//(num, u::T) where T<:EnergyUnit = T(num // val(u), -dim(u))
*(u::T, num) where T<:EnergyUnit = T(val(u) * num, dim(u))
/(u::T, num) where T<:EnergyUnit = T(val(u) / num, dim(u))
//(u::T, num) where T<:EnergyUnit = T(val(u) // num, dim(u))
^(u::T, num) where T<:EnergyUnit = T(val(u)^num, dim(u)*num)

==(u1::T, u2::T) where T<:EnergyUnit = val(u1) == val(u2) && dim(u1) == dim(u2)
==(u1::EnergyUnit, u2::EnergyUnit) = ==(promote(u1, u2)...)
isless(u1::T, u2::T) where T<:EnergyUnit = if __is_same_dimension(u1, u2)
    isless(val(u1), val(u2))
else
    (throw ∘ __diff_dimension_error)(u1, u2, "compare")
end
isless(u1::EnergyUnit, u2::EnergyUnit) = isless(promote(u1, u2)...)

sqrt(u::T) where T<:EnergyUnit = T(sqrt(val(u)), dim(u) // 2)
cbrt(u::T) where T<:EnergyUnit = T(cbrt(val(u)), dim(u) // 3)

dim(u::EnergyUnit) = u.dimension
val(u::EnergyUnit) = u.value
val(::Type{T}, u::EnergyUnit) where T<:EnergyUnit = convert(T, u).value
val(num::Number) = identity(num)
val(::Type{T}, num::Number) where T<:EnergyUnit = identity(num)
