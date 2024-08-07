# Copyright (c) 2024 Quan-feng WU <wuquanfeng@ihep.ac.cn>
# 
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

macro check_EU_dimension(expr, dim)
    return quote
        @assert EUdim($(esc(expr))) == $(esc(dim)) "Expected mass dimension of $(string($(esc(dim)))), but got $(EUdim($(esc(expr)))): $($(esc(expr)))."
    end
end

macro check_positive_value(expr)
    return quote
        @assert $(esc(expr)) > zero($(esc(expr))) "Expected positive value, but got $($(esc(expr)))."
    end
end

macro check_nonnegative_value(expr)
    return quote
        @assert $(esc(expr)) ≥ zero($(esc(expr))) "Expected non-negative value, but got $($(esc(expr)))."
    end
end
