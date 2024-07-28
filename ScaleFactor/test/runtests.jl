using NaturalUnits
using ScaleFactor

RD_test_a = 1e20
RD_result = scale_factor_to_time(Val(:RD), RD_test_a)
RD_rel_err = (
    time_to_scale_factor(Val(:RD), RD_result) - RD_test_a
) / RD_test_a
@show abs(RD_rel_err)

a_eq = 1/3403
a_ratio = 10

MR_result = scale_factor_to_time(Val(:MR), a_ratio * a_eq)
MR_rel_err = (
    time_to_scale_factor(Val(:MR), MR_result) - a_ratio * a_eq
) / (a_ratio * a_eq)
@show abs(MR_rel_err)

MR_result = scale_factor_to_time(Val(:MR), a_eq / a_ratio)
MR_rel_err = (
    time_to_scale_factor(Val(:MR), MR_result) - a_eq / a_ratio
) / (a_eq / a_ratio)
@show abs(MR_rel_err)
