__precompile__(true)

module DiffEqProblemLibrary

const file = @__FILE__

module ODEProblemLibrary
importodeproblems() = @eval ODEProblemLibrary include(joinpath(@__DIR__, "ode/ode_premade_problems.jl"))
end # module
module DAEProblemLibrary
importdaeproblems() = @eval DAEProblemLibrary include(joinpath(@__DIR__, "dae_premade_problems.jl"))
end # module
module DDEProblemLibrary
importddeproblems() = @eval DDEProblemLibrary include(joinpath(@__DIR__, "dde_premade_problems.jl"))
end # module
module SDEProblemLibrary
importsdeproblems() = @eval SDEProblemLibrary include(joinpath(@__DIR__, "sde_premade_problems.jl"))
end # module
module JumpProblemLibrary
importjumpproblems() = @eval JumpProblemLibrary include(joinpath(@__DIR__, "src/jump_premade_problems.jl"))
end # module

end # module
