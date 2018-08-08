module DiffEqProblemLibrary

module ODEProblemLibrary
function importodeproblems()
  @isdefined(prob_ode_linear) ||
  include(joinpath(@__DIR__, "ode/ode_premade_problems.jl"));
  nothing
end
end # module

module DAEProblemLibrary
function importdaeproblems()
  @isdefined(prob_dae_resrob) ||
  include(joinpath(@__DIR__, "dae_premade_problems.jl"));
  nothing
end
end # module

module DDEProblemLibrary
function importddeproblems()
  @isdefined(prob_dde_1delay) ||
  include(joinpath(@__DIR__, "dde_premade_problems.jl"));
  nothing
end
end # module

module SDEProblemLibrary
function importsdeproblems()
  @isdefined(prob_sde_wave) ||
  include(joinpath(@__DIR__, "sde_premade_problems.jl"));
  nothing
end
end # module

module JumpProblemLibrary
function importjumpproblems()
  @isdefined(prob_jump_dnarepressor) ||
  include(joinpath(@__DIR__, "jump_premade_problems.jl"));
  nothing
end
end # module

end # module
