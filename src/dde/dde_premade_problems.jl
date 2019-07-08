using DiffEqBase
# DDE Example Problems

# examples with constant delays
export
  # DDEs with 1 constant delay
  prob_dde_1delay_ip, prob_dde_1delay_oop, prob_dde_1delay_scalar,
  prob_dde_1delay_long_ip, prob_dde_1delay_long_oop, prob_dde_1delay_long_scalar,
  # DDEs with 2 constant delays
  prob_dde_2delays_ip, prob_dde_2delays_oop, prob_dde_2delays_scalar,
  prob_dde_2delays_long_ip, prob_dde_2delays_long_oop, prob_dde_2delays_long_scalar

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
  prob_dde_DDETST_F1, prob_dde_DDETST_F2, prob_dde_DDETST_F3, prob_dde_DDETST_F4, prob_dde_DDETST_F5,
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

# deprecations for problems with constant delays
for p in (:(1delay), :(1delay_long), :(2delays), :(2delays_long))
  @eval begin
    Base.@deprecate_binding $(Symbol(:prob_dde_, p)) $(Symbol(:prob_dde_constant_, p, :_ip))
    Base.@deprecate_binding $(Symbol(:prob_dde_, p, :_notinplace)) $(Symbol(:prob_dde_constant_, p, :_oop))
    Base.@deprecate_binding $(Symbol(:prob_dde_, p, :_scalar_notinplace)) $(Symbol(:prob_dde_constant_, p, :_scalar))

    Base.@deprecate $(Symbol(:build_prob_dde_, p))(u₀, t = u₀) remake_dde_constant_u0_tType($(Symbol(:prob_dde_constant_, p, :_ip)), [u₀], typeof(t))
    Base.@deprecate $(Symbol(:build_prob_dde_, p, :_notinplace))(u₀, t = u₀) remake_dde_constant_u0_tType($(Symbol(:prob_dde_constant_, p, :_oop)), [u₀], typeof(t))
    Base.@deprecate $(Symbol(:build_prob_dde_, p, :_scalar_notinplace))(u₀, t = u₀) remake_dde_constant_u0_tType($(Symbol(:prob_dde_constant_, p, :_scalar)), u₀, typeof(t))
  end
end

# deprecations for DDETST problems
Base.@deprecate_binding prob_dde_mackey prob_dde_DDETST_A1
Base.@deprecate_binding prob_dde_wheldon prob_dde_DDETST_A2
Base.@deprecate_binding prob_dde_neves1 prob_dde_DDETST_B1
Base.@deprecate_binding prob_dde_neves_thompson prob_dde_DDETST_B2
Base.@deprecate_binding prob_dde_paul1 prob_dde_DDETST_C1
Base.@deprecate_binding prob_dde_paul2 prob_dde_DDETST_C2
Base.@deprecate_binding prob_dde_mahaffy1 prob_dde_DDETST_C3
Base.@deprecate_binding prob_dde_mahaffy2 prob_dde_DDETST_C4
Base.@deprecate_binding prob_dde_neves2 prob_dde_DDETST_D1
Base.@deprecate_binding prob_dde_gatica prob_dde_DDETST_D2
