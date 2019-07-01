using DiffEqBase, ParameterizedFunctions, DiffEqOperators, Random, LinearAlgebra
using Markdown

Random.seed!(100)

#ODE Example Problems
export prob_ode_linear, prob_ode_bigfloatlinear, prob_ode_2Dlinear,
       prob_ode_large2Dlinear, prob_ode_bigfloat2Dlinear, prob_ode_rigidbody,
       prob_ode_2Dlinear_notinplace, prob_ode_vanderpol, prob_ode_vanderpol_stiff,
       prob_ode_lorenz, prob_ode_rober, prob_ode_threebody, prob_ode_mm_linear, prob_ode_pleiades,
       prob_ode_brusselator_1d, prob_ode_brusselator_2d, prob_ode_orego,
       prob_ode_hires, prob_ode_pollution, prob_ode_filament,
       prob_ode_nonlinchem,
       SolverDiffEq, filament_prob

include("ode_linear_prob.jl")
include("ode_simple_nonlinear_prob.jl")
include("brusselator_prob.jl")
include("pollution_prob.jl")
include("filament_prob.jl")
include("nonlinchem.jl")
