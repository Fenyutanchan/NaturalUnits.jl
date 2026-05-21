# Copyright (c) 2024 Quan-feng WU <wuquanfeng@ihep.ac.cn>
# 
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

macro check_EU_dimension(dim, exprs...)
    checks = map(exprs) do e
        quote
            local _x = $(esc(e))
            EUdim(_x) == _d || throw(
                DimensionMismatch("Expected mass dimension $(_d), got $(EUdim(_x)): $(_x).")
            )
        end
    end
    quote
        local _d = $(esc(dim))
        $(checks...)
    end
end

macro check_positive_value(exprs...)
    checks = map(exprs) do e
        quote
            local _x = $(esc(e))
            _x > zero(_x) || throw(
                ArgumentError("Expected positive value, got $(_x).")
            )
        end
    end
    Expr(:block, checks...)
end

macro check_nonnegative_value(exprs...)
    checks = map(exprs) do e
        quote
            local _x = $(esc(e))
            _x ≥ zero(_x) || throw(
                ArgumentError("Expected non-negative value, got $(_x).")
            )
        end
    end
    Expr(:block, checks...)
end

# Re-type the `value` field of an `EnergyUnit` without changing its prefix or
# dimension. Useful in generic code that needs to switch numeric precision
# (e.g. `Float64` → `BigFloat`) without knowing the concrete prefix of `u`.
convert_EnergyUnit_value_type(::Type{S}, u::U) where {S, U<:EnergyUnit} =
    _wrapper(U){S}(EUval(u), EUdim(u))
