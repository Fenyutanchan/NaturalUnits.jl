# Copyright (c) 2024 Quan-feng WU <wuquanfeng@ihep.ac.cn>
# 
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

# =============================================================================
# Abstract type
# =============================================================================

abstract type EnergyUnit{T<:Number} end

# =============================================================================
# Concrete-type generation
# =============================================================================

# Defines a concrete `EnergyUnit{T}` subtype `name` together with its
# constructors, `_wrapper`, same-prefix value-type conversion, and same-prefix
# promotion rule. Used both for the bare `eV` and for the prefixed siblings.
macro __define_EU_type(name)
    return quote
        export $(name)

        struct $(name){T<:Number} <: EnergyUnit{T}
            value::T
            dimension::Rational{Int}

            $(name){T}(value, dimension) where {T<:Number} =
                new{T}(convert(T, value), rationalize(Int, dimension))
        end

        $(name)() = $(name){Int}(1, 1)
        $(name)(value::T) where {T<:Number} = $(name){T}(value, 1)
        $(name)(value::T, dimension) where {T<:Number} =
            $(name){T}(value, dimension)

        # `_wrapper(U)` strips the value-type parameter from a concrete
        # `EnergyUnit` subtype, so we can reconstruct the same prefix with a
        # possibly different value type (e.g. after promotion in `num * u`).
        _wrapper(::Type{<:$(name)}) = $(name)

        # Same-prefix conversion that only changes the value type.
        convert(::Type{$(name){S}}, u::$(name)) where {S} =
            $(name){S}(EUval(u), EUdim(u))

        # Same-prefix promotion of value types.
        promote_rule(::Type{$(name){T}}, ::Type{$(name){S}}) where {T, S} =
            $(name){promote_type(T, S)}
    end |> esc
end

@__define_EU_type eV

# Defines a `XeV` prefixed sibling of `eV`, where `num` is the conversion
# factor `1 XeV = num eV`. Adds cross-prefix conversions to/from `eV` and
# attaches a short docstring to the new type.
macro generate_XeV(head, num)
    XeV = Symbol(head, "eV")
    doc = "`$(XeV)` — energy unit with `1 $(XeV) = $(num) eV`."
    return quote
        NaturalUnits.@__define_EU_type $(XeV)

        # Cross-prefix conversions to/from `eV` (untyped and typed targets).
        convert(::Type{eV}, u::$(XeV)) =
            eV(EUval(u) * $(num)^EUdim(u), EUdim(u))
        convert(::Type{eV{S}}, u::$(XeV)) where {S} =
            eV{S}(EUval(u) * $(num)^EUdim(u), EUdim(u))
        convert(::Type{$(XeV)}, u::eV) =
            $(XeV)(EUval(u) / $(num)^EUdim(u), EUdim(u))
        convert(::Type{$(XeV){S}}, u::eV) where {S} =
            $(XeV){S}(EUval(u) / $(num)^EUdim(u), EUdim(u))

        @doc $(doc) $(XeV)
    end |> esc
end

@generate_XeV k 1e3
@generate_XeV M 1e6
@generate_XeV G 1e9
@generate_XeV T 1e12

# =============================================================================
# Field accessors
# =============================================================================

EUdim(u::EnergyUnit) = u.dimension
EUdim(num::Number) = 0
EUval(u::EnergyUnit) = u.value
EUval(::Type{T}, u::EnergyUnit) where {T<:EnergyUnit} = convert(T, u).value
EUval(num::Number) = identity(num)
EUval(::Type{T}, num::Number) where {T<:EnergyUnit} = identity(num)
EUval(::Type{T}) where {T<:EnergyUnit} = x -> EUval(T, x)

# =============================================================================
# Generic conversion and promotion
# =============================================================================

convert(::Type{T}, u::T) where {T<:EnergyUnit} = identity(u)
convert(::Type{T}, num::Number) where {T<:EnergyUnit} = T(num, 0)
convert(::Type{EnergyUnit}, num::Number) = eV(num, 0)
convert(T::Type{<:EnergyUnit}, u::EnergyUnit) = convert(T, convert(eV, u))

# Same-prefix promotion (eV) and the cross-prefix fallback.
promote_rule(::Type{<:EnergyUnit{T}}, ::Type{<:EnergyUnit{S}}) where {T, S} =
    eV{promote_type(T, S)}
# Fallback for `UnionAll` arguments (e.g. `promote_rule(MeV, eV)` with no
# value-type parameter), preserved for backwards compatibility.
promote_rule(::Type{<:EnergyUnit}, ::Type{<:EnergyUnit}) = eV

# =============================================================================
# Identity elements
# =============================================================================

one(u::U) where {U<:EnergyUnit} = U(one(EUval(u)), 0)
one(::Type{T}) where {T<:EnergyUnit} = T(1, 0)

oneunit(u::U) where {U<:EnergyUnit} = U(one(EUval(u)), EUdim(u))
oneunit(::Type{T}) where {T<:EnergyUnit} = T(1, 1)

zero(u::U) where {U<:EnergyUnit} = U(zero(EUval(u)), EUdim(u))
# `zero(::Type{T})` is intentionally not defined: the additive identity depends
# on the dimension (a field, not a type parameter), so a single canonical zero
# does not exist for the bare type. Use `zero(u)` with an instance, or `iszero`
# to test for the additive identity at any dimension.

# =============================================================================
# Arithmetic
# =============================================================================

__is_same_dimension(u1::EnergyUnit, u2::EnergyUnit) = EUdim(u1) == EUdim(u2)
__diff_dimension_error(u1::EnergyUnit, u2::EnergyUnit, operate::String) = ArgumentError("Cannot $operate energy units with different dimensions: $(EUdim(u1)) and $(EUdim(u2)).")

# --- Addition and subtraction ---

+(u1::T, u2::T) where {T<:EnergyUnit} = if __is_same_dimension(u1, u2)
    T(EUval(u1) + EUval(u2), EUdim(u1))
else
    throw(__diff_dimension_error(u1, u2, "add"))
end
-(u1::T, u2::T) where {T<:EnergyUnit} = if __is_same_dimension(u1, u2)
    T(EUval(u1) - EUval(u2), EUdim(u1))
else
    throw(__diff_dimension_error(u1, u2, "subtract"))
end
+(u1::EnergyUnit, u2::EnergyUnit) = +(promote(u1, u2)...)
-(u1::EnergyUnit, u2::EnergyUnit) = -(promote(u1, u2)...)
-(u::U) where {U<:EnergyUnit} = U(-EUval(u), EUdim(u))

# Numeric zero acts as the additive identity at any dimension; non-zero numbers
# are converted to a dimensionless energy unit and will raise a dimension
# mismatch unless `u` is itself dimensionless.
+(u::U, x::Number) where {U<:EnergyUnit} = iszero(x) ? u : u + convert(U, x)
+(x::Number, u::U) where {U<:EnergyUnit} = iszero(x) ? u : convert(U, x) + u
-(u::U, x::Number) where {U<:EnergyUnit} = iszero(x) ? u : u - convert(U, x)
-(x::Number, u::U) where {U<:EnergyUnit} = iszero(x) ? -u : convert(U, x) - u

# --- Multiplication, division, exponentiation ---

*(u1::T, u2::T) where {T<:EnergyUnit} = T(EUval(u1) * EUval(u2), EUdim(u1) + EUdim(u2))
/(u1::T, u2::T) where {T<:EnergyUnit} = T(EUval(u1) / EUval(u2), EUdim(u1) - EUdim(u2))
//(u1::T, u2::T) where {T<:EnergyUnit} = T(EUval(u1) // EUval(u2), EUdim(u1) - EUdim(u2))

*(u1::EnergyUnit, u2::EnergyUnit) = *(promote(u1, u2)...)
/(u1::EnergyUnit, u2::EnergyUnit) = /(promote(u1, u2)...)
//(u1::EnergyUnit, u2::EnergyUnit) = //(promote(u1, u2)...)

# Mixed scalar / unit operators: use `_wrapper(U)` so the result's value type
# follows ordinary numeric promotion instead of being clamped to `U`'s
# parameter (which would silently truncate `Float * eV{Int}` etc.).
*(num::Number, u::U) where {U<:EnergyUnit} = _wrapper(U)(num * EUval(u), EUdim(u))
*(u::U, num::Number) where {U<:EnergyUnit} = _wrapper(U)(EUval(u) * num, EUdim(u))
/(num::Number, u::U) where {U<:EnergyUnit} = _wrapper(U)(num / EUval(u), -EUdim(u))
/(u::U, num::Number) where {U<:EnergyUnit} = _wrapper(U)(EUval(u) / num, EUdim(u))
//(num::Number, u::U) where {U<:EnergyUnit} = _wrapper(U)(num // EUval(u), -EUdim(u))
//(u::U, num::Number) where {U<:EnergyUnit} = _wrapper(U)(EUval(u) // num, EUdim(u))

^(u::U, num) where {U<:EnergyUnit} = _wrapper(U)(EUval(u)^num, EUdim(u) * num)
inv(u::U) where {U<:EnergyUnit} = _wrapper(U)(inv(EUval(u)), -EUdim(u))

# =============================================================================
# Comparison and hashing
# =============================================================================

==(u1::T, u2::T) where {T<:EnergyUnit} = EUval(u1) == EUval(u2) && EUdim(u1) == EUdim(u2)
==(u1::EnergyUnit, u2::EnergyUnit) = ==(promote(u1, u2)...)

isequal(u1::T, u2::T) where {T<:EnergyUnit} = isequal(EUval(u1), EUval(u2)) && EUdim(u1) == EUdim(u2)
isequal(u1::EnergyUnit, u2::EnergyUnit) = isequal(promote(u1, u2)...)

hash(u::EnergyUnit, h::UInt) = hash(EUdim(u), hash(EUval(eV, u), hash(:NaturalUnits_EnergyUnit, h)))

isless(u1::T, u2::T) where {T<:EnergyUnit} = if __is_same_dimension(u1, u2)
    isless(EUval(u1), EUval(u2))
else
    throw(__diff_dimension_error(u1, u2, "compare"))
end
isless(u1::EnergyUnit, u2::EnergyUnit) = isless(promote(u1, u2)...)

# =============================================================================
# Math functions
# =============================================================================

abs(u::U) where {U<:EnergyUnit} = _wrapper(U)(abs(EUval(u)), EUdim(u))
abs2(u::U) where {U<:EnergyUnit} = _wrapper(U)(abs2(EUval(u)), 2 * EUdim(u))
sqrt(u::U) where {U<:EnergyUnit} = _wrapper(U)(sqrt(EUval(u)), EUdim(u) // 2)
cbrt(u::U) where {U<:EnergyUnit} = _wrapper(U)(cbrt(EUval(u)), EUdim(u) // 3)
real(u::U) where {U<:EnergyUnit} = _wrapper(U)(real(EUval(u)), EUdim(u))
imag(u::U) where {U<:EnergyUnit} = _wrapper(U)(imag(EUval(u)), EUdim(u))
conj(u::U) where {U<:EnergyUnit} = _wrapper(U)(conj(EUval(u)), EUdim(u))
angle(u::EnergyUnit) = angle(EUval(u))

# =============================================================================
# Predicates
# =============================================================================

isinf(u::EnergyUnit) = isinf(EUval(u)) || isinf(EUdim(u))
isnan(u::EnergyUnit) = isnan(EUval(u)) || isnan(EUdim(u))
iszero(u::EnergyUnit) = iszero(EUval(u))

# =============================================================================
# Broadcasting (treat a single unit as a scalar, like `Number`)
# =============================================================================

broadcastable(u::EnergyUnit) = Ref(u)
