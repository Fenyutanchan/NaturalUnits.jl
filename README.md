# NaturalUnits.jl

A lightweight Julia package for working with **natural units** ($\hbar = c = k_B = 1$) in high-energy physics. Physical quantities are expressed in units of electron-volts (eV) with automatic dimensional tracking through all arithmetic operations.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/Fenyutanchan/NaturalUnits.jl.git")
```

Or for development:

```julia
Pkg.develop(url="https://github.com/Fenyutanchan/NaturalUnits.jl.git")
```

## Quick Start

```julia
using NaturalUnits

# Create energy quantities with automatic dimensional tracking
E = GeV(125.0)          # 125 GeV, dimension = 1
m = GeV(125.0)          # mass, dimension = 1
λ = GeV(1.0, -1)        # length, dimension = -1 (GeV⁻¹)

# Arithmetic propagates dimensions
E² = E * E              # GeV(15625.0, 2)
E_ratio = E / m         # GeV(1.0, 0) — dimensionless, but still an EnergyUnit
λ² = λ^2                # GeV(1.0, -2)

# Mixed-prefix operations promote through `eV`
keV(1) + MeV(1)         # eV(1.001e6, 1)
```

## Type System

All energy units are parameterized by the value type:

```julia
eV{T<:Number} <: EnergyUnit{T}
keV{T<:Number} <: EnergyUnit{T}
MeV{T<:Number} <: EnergyUnit{T}
GeV{T<:Number} <: EnergyUnit{T}
TeV{T<:Number} <: EnergyUnit{T}
```

The constructors infer `T` from the value:

```julia
eV(3)             # eV{Int}(3, 1//1)
eV(3.0)           # eV{Float64}(3.0, 1//1)
eV(3.0, 2)        # eV{Float64}(3.0, 2//1)
```

The `dimension` field is stored as `Rational{Int}`, so fractional dimensions arising from `sqrt`, `cbrt`, etc. are exact:

```julia
sqrt(GeV(4.0, 2))   # GeV(2.0, 1//1)
cbrt(GeV(8.0, 3))   # GeV(2.0, 1//1)
sqrt(GeV(1.0))      # GeV(1.0, 1//2)
```

Concrete unit types are `isbits`, which means no boxing in tight loops:

```julia
isbitstype(GeV{Float64})   # true
```

## SI ↔ Natural Unit Conversions

`NaturalUnit{T}` provides a property-based interface for converting 1 SI unit into natural units:

```julia
nu = NaturalUnit(GeV)

nu.m      # 1 meter  in GeV⁻¹  ≈ 5.068e15
nu.cm     # 1 cm     in GeV⁻¹
nu.s      # 1 second in GeV⁻¹  ≈ 1.519e24
nu.kg     # 1 kg     in GeV    ≈ 5.610e26
nu.g      # 1 g      in GeV
nu.J      # 1 Joule  in GeV    ≈ 6.242e9
nu.K      # 1 Kelvin in GeV    ≈ 8.617e-14

nu.G_N    # Newton's gravitational constant in GeV⁻²
nu.M_Pl   # Reduced Planck mass ≈ 2.435e18 GeV
nu.m_Pl   # Planck mass         ≈ 1.221e19 GeV
```

### Typical usage: SI → natural

```julia
nu = NaturalUnit(GeV)

# Express a cross section (units of area) in GeV⁻²
σ = 1e-38 * nu.cm^2     # σ in GeV⁻²

# Express a lifetime in GeV⁻¹
τ = 1e-12 * nu.s        # τ in GeV⁻¹

# Express a temperature in GeV
T = 2.725 * nu.K        # CMB temperature in GeV
```

### Typical usage: natural → SI

```julia
nu = NaturalUnit(GeV)

# A decay width of 1 GeV corresponds to a lifetime of:
τ = 1.0 / GeV(1.0)              # in GeV⁻¹
τ_seconds = EUval(τ) / EUval(nu.s)   # back to seconds
```

## API Reference

### Energy unit types

| Type | Description |
|------|-------------|
| `EnergyUnit{T}` | Abstract supertype for all energy unit types |
| `eV{T}` | Electron-volt |
| `keV{T}` | Kiloelectron-volt ($10^{3}$ eV) |
| `MeV{T}` | Megaelectron-volt ($10^{6}$ eV) |
| `GeV{T}` | Gigaelectron-volt ($10^{9}$ eV) |
| `TeV{T}` | Teraelectron-volt ($10^{12}$ eV) |

### Constructors

```julia
eV()              # eV{Int}(1, 1)
eV(val)           # eV{typeof(val)}(val, 1)
eV(val, dim)      # eV{typeof(val)}(val, rationalize(dim))
```

Same pattern for `keV`, `MeV`, `GeV`, `TeV`.

### Accessors

| Function | Description |
|----------|-------------|
| `EUdim(u)` | Mass dimension of `u` as `Rational{Int}` (returns `0` for plain `Number`) |
| `EUval(u)` | Numeric value of `u` |
| `EUval(T, u)` | Value of `u` after conversion to unit type `T` |
| `EUval(T)` | Curried form: returns `x -> EUval(T, x)` |

### Dimensional arithmetic

Arithmetic operators automatically track mass dimensions:

| Operation | Dimension rule |
|-----------|----------------|
| `a + b`, `a - b` | Dimensions must match (throws on mismatch) |
| `a * b` | Dimensions add |
| `a / b`, `a // b` | Dimensions subtract |
| `a^n` | Dimension multiplied by `n` |
| `inv(a)` | Dimension negated |
| `sqrt(a)` | Dimension halved |
| `cbrt(a)` | Dimension divided by 3 |

Mixed-prefix operations promote through `eV{promote_type(T, S)}`:

```julia
GeV(1) + MeV(1000.0)    # eV(2.0e9, 1//1)
```

Mixed scalar / unit operations follow ordinary numeric promotion on the value
type:

```julia
2.5 * eV{Int}(3, 1)     # eV{Float64}(7.5, 1//1)
```

### Conversion

```julia
convert(GeV, MeV(1000.0))                          # GeV(1.0, 1//1)
convert(GeV{Float32}, MeV(1000.0))                 # GeV{Float32}(1.0f0, 1//1)
convert_EnergyUnit_value_type(BigFloat, GeV(5.0))  # GeV{BigFloat}(5.0, 1//1)
```

- `convert(T, u)` may change the **prefix** (rescales the value).
- `convert_EnergyUnit_value_type(S, u)` only changes the **value type**, keeping
  the prefix and dimension unchanged. Useful in generic code that wants to switch
  numeric precision without knowing the concrete prefix.

### Comparison and hashing

`==`, `isequal`, `isless`, `hash` are defined. Equality and ordering promote
across prefixes (so `keV(1) == eV(1000)` is `true`), and `isless` errors on
dimension mismatch.

### Mathematical functions

`abs`, `abs2`, `sqrt`, `cbrt`, `real`, `imag`, `conj`, `angle`, `isinf`,
`isnan`, `iszero` are all supported and propagate dimensions appropriately.

### Broadcasting

`EnergyUnit` values behave as **scalars** in broadcasting, like `Number`:

```julia
u   = eV(3.0, 1)
arr = [1.0, 2.0, 3.0]

EUval.(u)     # 3.0
arr .* u      # [eV(3.0, 1), eV(6.0, 1), eV(9.0, 1)]
```

(The unit itself is not iterable; `for x in u` will not work, as is the case
for `Number`.)

### Defining new prefixed units

The `@generate_XeV` macro defines a new prefix sibling of `eV` and attaches
conversions to/from `eV`. The package itself uses it to define `keV`, `MeV`,
`GeV`, `TeV`:

```julia
@generate_XeV P 1e15    # defines PeV with 1 PeV = 1e15 eV
PeV(1.0) + TeV(1000.0)  # promotes to eV
```

### Custom NaturalUnit properties

Extend `NaturalUnit` with your own conversion functions:

```julia
function __femtometer(nu::NaturalUnit)
    return nu.m * 1e-15
end
add_property_function(:fm, __femtometer)

nu = NaturalUnit(GeV)
nu.fm    # 1 fm in GeV⁻¹
```

### Assertion macros

All three macros are **variadic** and evaluate each argument exactly once
(so they are safe to call on expressions with side effects). On failure they
throw a real exception (not an `@assert`-style `AssertionError` that can be
elided by the compiler).

| Macro | Signature | Throws |
|-------|-----------|--------|
| `@check_EU_dimension dim exprs...` | Assert every `EUdim(expr) == dim` | `DimensionMismatch` |
| `@check_positive_value exprs...` | Assert every `expr > zero(expr)` | `ArgumentError` |
| `@check_nonnegative_value exprs...` | Assert every `expr ≥ zero(expr)` | `ArgumentError` |

Examples:

```julia
@check_EU_dimension 2 u v w        # all three must have dimension 2
@check_positive_value E1 E2        # both must be positive
```

## Physical Constants

The package uses CODATA 2019 exact SI values:

| Constant | Value |
|----------|-------|
| $h$ | $6.62607015 \times 10^{-34}$ J·s |
| $\hbar$ | $h / (2\pi)$ |
| $c$ | $299\,792\,458$ m/s |
| $k_B$ | $1.380649 \times 10^{-23}$ J/K |
| $e$ | $1.602176634 \times 10^{-19}$ C |

## Testing

```julia
using Pkg
Pkg.test("NaturalUnits")
```

## License

MIT License.
Copyright (c) 2024 Quan-feng WU.
