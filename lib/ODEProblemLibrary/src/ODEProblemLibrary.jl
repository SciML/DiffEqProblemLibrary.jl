module ODEProblemLibrary

using DiffEqBase
using Latexify
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

using LinearAlgebra
using Markdown
using Random

import RuntimeGeneratedFunctions
RuntimeGeneratedFunctions.init(@__MODULE__)

Random.seed!(100)

#ODE Example Problems
export prob_ode_linear, prob_ode_bigfloatlinear, prob_ode_2Dlinear,
       prob_ode_large2Dlinear, prob_ode_bigfloat2Dlinear, prob_ode_rigidbody,
       prob_ode_2Dlinear_notinplace, prob_ode_vanderpol, prob_ode_vanderpol_stiff,
       prob_ode_lotkavolterra, prob_ode_fitzhughnagumo,
       prob_ode_rober, prob_ode_threebody, prob_ode_mm_linear, prob_ode_pleiades,
       prob_ode_brusselator_1d, prob_ode_brusselator_2d, prob_ode_orego,
       prob_ode_hires, prob_ode_pollution, prob_ode_filament,
       prob_ode_nonlinchem,
       SolverDiffEq, filament_prob,
       thomas, lorenz, aizawa, dadras, chen, rossler, rabinovich_fabrikant, sprott,
       hindmarsh_rose,
       prob_ode_thomas, prob_ode_lorenz, prob_ode_aizawa, prob_ode_dadras,
       prob_ode_chen, prob_ode_rossler, prob_ode_rabinovich_fabrikant, prob_ode_sprott,
       prob_ode_hindmarsh_rose

include("ode_linear_prob.jl")
include("ode_simple_nonlinear_prob.jl")
include("brusselator_prob.jl")
include("pollution_prob.jl")
include("filament_prob.jl")
include("nonlinchem.jl")
include("strange_attractors.jl")

end # module
