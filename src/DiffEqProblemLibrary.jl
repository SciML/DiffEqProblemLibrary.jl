__precompile__(true)

module DiffEqProblemLibrary

module ODEProblemLibrary
importodeproblems() = include(joinpath(@__DIR__, "ode/ode_premade_problems.jl"))
end # module

module DAEProblemLibrary
importdaeproblems() = include(joinpath(@__DIR__, "dae_premade_problems.jl"))
end # module

module DDEProblemLibrary
importddeproblems() = include(joinpath(@__DIR__, "dde_premade_problems.jl"))
end # module

module SDEProblemLibrary
importsdeproblems() = include(joinpath(@__DIR__, "sde_premade_problems.jl"))
end # module

module JumpProblemLibrary
importjumpproblems() = include(joinpath(@__DIR__, "jump_premade_problems.jl"))
end # module

end # module
