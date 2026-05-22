# Copyright (c) 2026 Quan-feng WU <wuquanfeng@ihep.ac.cn>
# 
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

import Base: log, log2, log10, exp, exp2, exp10

export unwrap_dimensionless_EU

function unwrap_dimensionless_EU(u::EnergyUnit{T})::Union{T, EnergyUnit{T}} where {T<:Number}
    return iszero(EUdim(u)) ? EUval(u) : u
end
unwrap_dimensionless_EU(u::Number) = identity(u)

function _unwrap_dimensionless_EU(u::EnergyUnit)
    iszero(EUdim(u)) || throw(ArgumentError("Energy unit must be dimensionless!"))
    return EUval(u)
end

log(u::EnergyUnit) = _unwrap_dimensionless_EU(u) |> log
log2(u::EnergyUnit) = _unwrap_dimensionless_EU(u) |> log2
log10(u::EnergyUnit) = _unwrap_dimensionless_EU(u) |> log10
exp(u::EnergyUnit) = _unwrap_dimensionless_EU(u) |> exp
exp2(u::EnergyUnit) = _unwrap_dimensionless_EU(u) |> exp2
exp10(u::EnergyUnit) = _unwrap_dimensionless_EU(u) |> exp10
