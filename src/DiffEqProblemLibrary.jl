__precompile__(true)

module DiffEqProblemLibrary

module ODEProblemLibrary
importodeproblems() =
@isdefined(prob_ode_linear) || include(joinpath(@__DIR__, "ode/ode_premade_problems.jl"))
end # module

module DAEProblemLibrary
importdaeproblems() =
@isdefined(prob_dae_resrob) || include(joinpath(@__DIR__, "dae_premade_problems.jl"))
end # module

module DDEProblemLibrary
importddeproblems() =
@isdefined(prob_dde_1delay) || include(joinpath(@__DIR__, "dde_premade_problems.jl"))
end # module

module SDEProblemLibrary
importsdeproblems() =
@isdefined(prob_sde_wave) || include(joinpath(@__DIR__, "sde_premade_problems.jl"))
end # module

module JumpProblemLibrary
importjumpproblems() =
@isdefined(prob_jump_dnarepressor) || include(joinpath(@__DIR__, "jump_premade_problems.jl"))
end # module

end # module
