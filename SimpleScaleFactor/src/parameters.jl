include("parameter_constructors.jl")

const __parameter_constructor_dict = OrderedDict{String, Function}(
    "H_0" => __set_H_0!,

    "a_mΛ" => __set_a_mΛ!,
    "a_eq" => __set_a_eq!,

    "t_0" => __set_t_0!,
    "t_mΛ" => __set_t_mΛ!,

    "η_star" => __set_η_star!,
    "η_eq" => __set_η_eq!,

    "η_m" => __set_η_m!,
    "t_m" => __set_t_m!,

    "η_r" => __set_η_r!,
    "t_r" => __set_t_r!,

    "t_EW" => __set_t_EW!
)
const __required_defaults = Dict{String, Any}(
    "a_0" => 1,
    "tolerance" => 1e-6,
    "Ω_m" => .3153,
    "Ω_r" => 9.02e-5,
    "Ω_Λ" => .6847,
    "h" => .674,
    "t_eq" => 51100 * (365 + 1/4) * 24 * 60^2 * NaturalUnit(MeV).s,
    "T_EW" => GeV(1e2),
    "g_star" => 106.75,
    "EU" => MeV,
    "NU" => NaturalUnit(MeV)
)
const __forbidden_keys = keys(__parameter_constructor_dict)
const __required_keys = keys(__required_defaults)
const __all_keys = union(__forbidden_keys, __required_keys)

set_parameters() = set_parameters(__required_defaults)
function set_parameters(input_param_dict::Dict)
    forbidden_keys = intersect(__forbidden_keys, keys(input_param_dict))
    isempty(forbidden_keys) ||
        @warn "Forbidden keys (will be ignored): $(join(forbidden_keys, ", "))"

    other_keys = setdiff(keys(input_param_dict), __all_keys)
    isempty(other_keys) ||
        @warn "Unsupported keys (will be ignored): $(join(other_keys, ", "))"

    all_parameters = Dict{String, Any}()
    missing_keys = setdiff(__required_keys, keys(input_param_dict))
    for key ∈ missing_keys
        value = __required_defaults[key]
        @warn "Missing required key: $key. Using default value: $value."
        all_parameters[key] = value
    end
    for key ∈ setdiff(__required_keys, missing_keys)
        all_parameters[key] = input_param_dict[key]
    end
    for key ∈ __forbidden_keys
        __parameter_constructor_dict[key](all_parameters)
    end

    return all_parameters
end
