module DDEProblemLibrary

using DiffEqBase
# DDE Example Problems

# examples with constant delays
export
# DDEs with 1 constant delay
      prob_dde_constant_1delay_ip, prob_dde_constant_1delay_oop,
      prob_dde_constant_1delay_scalar,
      prob_dde_constant_1delay_long_ip, prob_dde_constant_1delay_long_oop,
      prob_dde_constant_1delay_long_scalar,
# DDEs with 2 constant delays
      prob_dde_constant_2delays_ip, prob_dde_constant_2delays_oop,
      prob_dde_constant_2delays_scalar,
      prob_dde_constant_2delays_long_ip, prob_dde_constant_2delays_long_oop,
      prob_dde_constant_2delays_long_scalar

# DDETST problems
export
# DDEs with time dependent delays
      prob_dde_DDETST_A1, prob_dde_DDETST_A2,
# DDEs with vanishing time dependent delays
      prob_dde_DDETST_B1, prob_dde_DDETST_B2,
# DDEs with state dependent delays
      prob_dde_DDETST_C1, prob_dde_DDETST_C2, prob_dde_DDETST_C3, prob_dde_DDETST_C4,
# DDEs with vanishing state dependent delays
      prob_dde_DDETST_D1, prob_dde_DDETST_D2,
# neutral DDEs with time dependent delays
      prob_dde_DDETST_E1, prob_dde_DDETST_E2,
# neutral DDEs with vanishing time dependent delays
      prob_dde_DDETST_F1, prob_dde_DDETST_F2, prob_dde_DDETST_F3, prob_dde_DDETST_F4,
      prob_dde_DDETST_F5,
# neutral DDEs with state dependent delays
      prob_dde_DDETST_G1, prob_dde_DDETST_G2,
# neutral DDEs with vanishing state dependent delays
      prob_dde_DDETST_H1, prob_dde_DDETST_H2, prob_dde_DDETST_H3, prob_dde_DDETST_H4

# RADAR5 problems
export prob_dde_RADAR5_oregonator, prob_dde_RADAR5_robertson, prob_dde_RADAR5_waltman

# QS exaple
export prob_dde_qs

include("constant_delays.jl")
include("ddetst.jl")
include("qs.jl")
include("radar5.jl")

end # module
