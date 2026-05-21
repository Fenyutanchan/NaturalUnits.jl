# Copyright (c) 2024 Quan-feng WU <wuquanfeng@ihep.ac.cn>
# 
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

module SIConstants
    # CODATA / SI defining constants in base SI units.
    const h   = 6.62607015e-34       # J·s   (Planck constant)
    const ħ   = h / (2 * π)          # J·s   (reduced Planck constant)
    const c   = 299792458            # m/s   (speed of light in vacuum)
    const k_B = 1.380649e-23         # J/K   (Boltzmann constant)
    const e   = 1.602176634e-19      # C     (elementary charge)
    const G_N = 6.67430e-11          # m^3·kg^-1·s^-2  (Newtonian constant of gravitation, CODATA 2018)
end