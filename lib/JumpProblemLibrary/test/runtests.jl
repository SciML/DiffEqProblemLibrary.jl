using Pkg
using JumpProblemLibrary
using SafeTestsets, Test

const TEST_GROUP = get(ENV, "GROUP", "All")

function activate_qa_env()
    Pkg.activate(joinpath(@__DIR__, "qa"))
    Pkg.develop(Pkg.PackageSpec(path = dirname(@__DIR__)))
    return Pkg.instantiate()
end

if TEST_GROUP == "Core" || TEST_GROUP == "All"
    @time @testset "Load Tests" begin
        @test JumpProblemLibrary isa Module
    end
    @time @safetestset "Network data" include("network_tests.jl")
end

# Quality assurance: no undefined exports, stale dependencies, etc.
if TEST_GROUP == "QA" || TEST_GROUP == "All"
    activate_qa_env()
    @time @safetestset "Aqua" include("qa/qa.jl")
end
