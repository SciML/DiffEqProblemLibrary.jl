__precompile__(false)

module DiffEqProblemLibrary

using DiffEqBase, ParameterizedFunctions, DiffEqOperators, DiffEqBiological

include("ode/ode_premade_problems.jl")
include("dae_premade_problems.jl")
include("dde_premade_problems.jl")
include("sde_premade_problems.jl")
include("jump_premade_problems.jl")

#ODE Example Problems
export prob_ode_linear, prob_ode_bigfloatlinear, prob_ode_2Dlinear,
       prob_ode_large2Dlinear, prob_ode_bigfloat2Dlinear, prob_ode_rigidbody,
       prob_ode_2Dlinear_notinplace, prob_ode_vanderpol, prob_ode_vanderpol_stiff,
       prob_ode_lorenz, prob_ode_rober, prob_ode_threebody, prob_ode_mm_linear, prob_ode_pleiades,
       prob_ode_brusselator_1d, prob_ode_brusselator_2d, prob_ode_orego,
       prob_ode_hires, prob_ode_pollution, prob_ode_filament,
       SolverDiffEq, filament_prob

 #SDE Example Problems
 export prob_sde_wave, prob_sde_linear, prob_sde_cubic, prob_sde_2Dlinear,
        prob_sde_lorenz, prob_sde_2Dlinear, prob_sde_additive,
        prob_sde_additivesystem, oval2ModelExample, prob_sde_bistable,
        prob_sde_bruss, prob_sde_oscilreact

 #SDE Stratonovich Example Problems
 export prob_sde_linear_stratonovich, prob_sde_2Dlinear_stratonovich

 #DAE Example Problems
 export prob_dae_resrob

# DDE Example Problems
# examples with constant delays
export prob_dde_1delay, prob_dde_1delay_notinplace, prob_dde_1delay_scalar_notinplace,
       prob_dde_2delays, prob_dde_2delays_notinplace, prob_dde_2delays_scalar_notinplace,
       prob_dde_1delay_long, prob_dde_1delay_long_notinplace,
       prob_dde_1delay_long_scalar_notinplace, prob_dde_2delays_long,
       prob_dde_2delays_long_notinplace, prob_dde_2delays_long_scalar_notinplace,
       prob_dde_mackey, prob_dde_wheldon, prob_dde_qs,
# examples with vanishing time dependent delays
       prob_ddde_neves1, prob_dde_neves_thompson,
# examples with state dependent delays
       prob_dde_paul1, prob_dde_paul2, prob_dde_mahaffy1, prob_dde_mahaffy2,
# examples with vanishing state dependent delays
       prob_neves2, prob_dde_gatica

# Jump Example Problems
export prob_jump_dnarepressor, prob_jump_constproduct, prob_jump_nonlinrxs,
# examples mixing mass action and constant rate jumps
       prob_jump_osc_mixed_jumptypes,
# examples used in published benchmarks / comparisions
       prob_jump_multistate, prob_jump_twentygenes, prob_jump_dnadimer_repressor,
# examples approximating diffusion by continuous time random walks
       prob_jump_diffnetwork

end # module
